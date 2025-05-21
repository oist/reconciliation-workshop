/*
    Copyright (C) 2017 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "corax/corax.h"

#define CORAX_PHYLIP_SEQUENTIAL 1
#define CORAX_PHYLIP_INTERLEAVED 2

static int
dfa_parse(corax_phylip_t *fd, corax_msa_t *msa, char *p, int seqno, int offset)
{
  int  j = 0;
  char c, m;

  char *seqdata = msa->sequence[seqno] + offset;

  /* read sequence data */
  while ((c = *p++))
  {
    m = (char)fd->chrstatus[(int)c];
    switch (m)
    {
    case 0:
      /* characters to be stripped */
      fd->stripped_count++;
      fd->stripped[(int)c]++;
      break;

    case 1:
      /* legal character */
      if (offset + j >= msa->length)
      {
        corax_set_error(CORAX_ERROR_PHYLIP_LONGSEQ,
                        "Sequence %d (%.100s) longer than expected",
                        seqno + 1,
                        msa->label[seqno]);
        return -1;
      }
      seqdata[j++] = c;
      break;

    case 2:
      /* fatal character */
      if (c >= 32)
      {
        corax_set_error(CORAX_ERROR_PHYLIP_ILLEGALCHAR,
                        "illegal character '%c' "
                        "on line %ld in the fasta file",
                        c,
                        fd->lineno);
      }
      else
      {
        corax_set_error(CORAX_ERROR_PHYLIP_UNPRINTABLECHAR,
                        "illegal unprintable character "
                        "%#.2x (hexadecimal) on line %ld "
                        "in the fasta file",
                        c,
                        fd->lineno);
      }
      return -1;

    case 3:
      /* silently stripped chars */
      break;
    }
  }
  return j;
}

static const int delimiters[4] = {' ', '\t', '\r', '\n'};
static int       get_headerlen(char *s)
{
  long min_len = strlen(s);
  int  i;

  for (i = 0; i < 4; ++i)
  {
    char *r = strchr(s, delimiters[i]);
    if (r && (r - s) < min_len) min_len = r - s;
  }

  return (int)min_len;
}

static char *reallocline(corax_phylip_t *fd, size_t newmaxsize)
{
  char *temp = (char *)malloc((size_t)newmaxsize * sizeof(char));
  if (!temp)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return NULL;
  }

  memcpy(temp, fd->line, fd->line_size * sizeof(char));
  free(fd->line);
  fd->line         = temp;
  fd->line_maxsize = newmaxsize;

  return temp;
}

static char *getnextline(corax_phylip_t *fd)
{
  size_t len = 0;

  fd->line_size = 0;

  /* read from file until newline or eof */
  while (fgets(fd->buffer, CORAX_LINEALLOC, fd->fp))
  {
    len = strlen(fd->buffer);

    if (fd->line_size + len > fd->line_maxsize)
      if (!reallocline(fd, fd->line_maxsize + CORAX_LINEALLOC)) return NULL;

    memcpy(fd->line + fd->line_size, fd->buffer, len * sizeof(char));
    fd->line_size += len;

    if (fd->buffer[len - 1] == '\n')
    {
#if 0
      if (line_size+1 > line_maxsize)
        reallocline(line_maxsize+1);

      line[line_size] = 0;
#else
      fd->line[fd->line_size - 1] = 0;
#endif

      return fd->line;
    }
  }

  if (!fd->line_size)
  {
    free(fd->line);
    fd->line = NULL;
    return NULL;
  }

  if (fd->line_size == fd->line_maxsize)
    if (!reallocline(fd, fd->line_maxsize + 1)) return NULL;

  fd->line[fd->line_size] = 0;
  return fd->line;
}

static int args_getint(const char *arg, int *len)
{
  int temp;
  *len = 0;

  int ret = sscanf(arg, "%d%n", &temp, len);
  if ((ret == 0) || (!*len)) return 0;

  return temp;
}

static int whitespace(char c)
{
  if (c == ' ' || c == '\t' || c == '\n' || c == '\r') return 1;
  return 0;
}

static int
parse_header(const char *line, int *seq_count, int *seq_len, int format)
{
  int len;

  /* read number of sequences */
  if (!(*seq_count = args_getint(line, &len)))
  {
    corax_set_error(CORAX_ERROR_PHYLIP_SYNTAX,
                    "Invalid number of sequences in header");
    return CORAX_FAILURE;
  }

  line += len;

  /* read sequence length */
  if (!(*seq_len = args_getint(line, &len)))
  {
    corax_set_error(CORAX_ERROR_PHYLIP_SYNTAX,
                    "Invalid sequence length in header");
    return CORAX_FAILURE;
  }

  line += len;

  /* go through all white spaces */
  while (*line && whitespace(*line)) ++line;

  /* if end of line then return successfully */
  if (!*line) return 1;

  /* otherwise, continue only if interleaved format specified, otherwise die */
  if (format == CORAX_PHYLIP_SEQUENTIAL) return 0;

  if (*line != 's' && *line != 'S' && *line != 'i' && *line != 'I') return 0;

  /* go through all white spaces */
  while (*line && whitespace(*line)) ++line;

  /* if end of line then return successfully */
  if (!*line) return 1;

  return 0;
}

static char *parse_oneline_sequence(corax_phylip_t *fd,
                                    corax_msa_t *   msa,
                                    char *          p,
                                    int             seqno,
                                    int             offset,
                                    int *           aln_len,
                                    int *           error)
{
  int j = 0;

  while (p && !j)
  {
    /* read data */
    if ((j = dfa_parse(fd, msa, p, seqno, offset)) == -1)
    {
      *error = 1;
      return NULL;
    }

    if (j)
    {
      if (!(*aln_len)) { *aln_len = j; }
      else if (*aln_len != j)
      {
        *error = 1;
        corax_set_error(CORAX_ERROR_PHYLIP_NONALIGNED,
                        "Sequence %d (%.100s) data out of alignment",
                        seqno + 1,
                        msa->label[seqno]);
        return NULL;
      }
    }
    else
      p = getnextline(fd);
  }

  return p;
}

CORAX_EXPORT corax_phylip_t *corax_phylip_open(const char *        filename,
                                               const unsigned int *map)
{
  int i;

  corax_phylip_t *fd = (corax_phylip_t *)malloc(sizeof(corax_phylip_t));
  if (!fd)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return NULL;
  }

  /* allocate space */
  fd->line         = NULL;
  fd->line_size    = 0;
  fd->line_maxsize = 0;

  fd->lineno = 0;

  fd->no = -1;

  fd->chrstatus = map;

  /* open file */
  fd->fp = fopen(filename, "r");
  if (!(fd->fp))
  {
    corax_set_error(
        CORAX_ERROR_FILE_OPEN, "Unable to open file (%s)", filename);
    free(fd);
    return NULL;
  }

  /* get filesize */
  if (fseek(fd->fp, 0, SEEK_END))
  {
    corax_set_error(
        CORAX_ERROR_FILE_SEEK, "Unable to seek in file (%s)", filename);
    fclose(fd->fp);
    free(fd);
    return NULL;
  }
  fd->filesize = ftell(fd->fp);

  rewind(fd->fp);

  /* reset stripped char frequencies */
  fd->stripped_count = 0;
  for (i = 0; i < 256; i++) fd->stripped[i] = 0;

  /* cache line */
  if (!getnextline(fd))
  {
    if (fd->line) free(fd->line);
    fclose(fd->fp);
    free(fd);
    return NULL;
  }

  fd->lineno = 1;

  return fd;
}

CORAX_EXPORT int corax_phylip_rewind(corax_phylip_t *fd)
{
  int i;

  rewind(fd->fp);

  /* reset stripped char frequencies */
  fd->stripped_count = 0;
  for (i = 0; i < 256; i++) fd->stripped[i] = 0;

  if (!getnextline(fd))
  {
    corax_set_error(CORAX_ERROR_FILE_SEEK, "Unable to rewind and cache data");
    return CORAX_FAILURE;
  }
  fd->lineno = 1;
  fd->no     = -1;

  return CORAX_SUCCESS;
}

CORAX_EXPORT void corax_phylip_close(corax_phylip_t *fd)
{
  fclose(fd->fp);
  if (fd->line) free(fd->line);
  free(fd);
}

CORAX_EXPORT corax_msa_t *corax_phylip_parse_interleaved(corax_phylip_t *fd)
{
  int  i;
  int  aln_len;
  int  sumlen;
  int  seqno;
  long headerlen;

  corax_msa_t *msa = (corax_msa_t *)malloc(sizeof(corax_msa_t));
  if (!msa)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return NULL;
  }

  /* read header */
  if (!parse_header(
          fd->line, &(msa->count), &(msa->length), CORAX_PHYLIP_INTERLEAVED))
  {
    free(msa);
    return NULL;
  }

  /* allocate msa placeholders */
  msa->sequence = (char **)calloc((size_t)(msa->count), sizeof(char *));
  msa->label    = (char **)calloc((size_t)(msa->count), sizeof(char *));
  if (!msa->label || !msa->sequence)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    corax_msa_destroy(msa);
    return NULL;
  }

  /* allocate sequence data placeholders */
  for (i = 0; i < msa->count; ++i)
  {
    msa->sequence[i] = (char *)malloc((size_t)(msa->length + 1) * sizeof(char));
    if (!msa->sequence[i])
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Unable to allocate enough memory.");
      corax_msa_destroy(msa);
      return NULL;
    }

    msa->sequence[i][msa->length] = 0;
  }

  /* assume the worst: whole sequence in one line + long label. *
   * Doing reallocate+copy really hurts on long alignments,     *
   * so we better waste some memory here to avoid reallocation. */
  reallocline(fd, (size_t)msa->length + 300);

  /* read sequences with headers */
  seqno     = 0;
  aln_len   = 0;
  int error = 0;
  while (1)
  {
    /* get next line */
    char *p = getnextline(fd);

    /* if no more lines break */
    if (!p) break;

    /* skip whitespace before sequence header */
    while (*p && whitespace(*p)) ++p;

    /* restart loop if blank line */
    if (!*p) continue;

    /* error if there are more sequences than specified */
    if (seqno == msa->count)
    {
      corax_set_error(CORAX_ERROR_PHYLIP_SYNTAX,
                      "Found at least %d sequences but expected %d",
                      seqno + 1,
                      msa->count);
      corax_msa_destroy(msa);
      return NULL;
    }

    /* find first delimiter after header */
    headerlen = get_headerlen(p);

    /* headerlen cannot be zero */
    assert(headerlen > 0);

    /* store sequence header */
    msa->label[seqno] = (char *)malloc((size_t)(headerlen + 1) * sizeof(char));
    if (!msa->label[seqno])
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Unable to allocate enough memory.");
      corax_msa_destroy(msa);
      return NULL;
    }
    memcpy(msa->label[seqno], p, (size_t)headerlen);
    msa->label[seqno][headerlen] = 0;

    p += headerlen;

    /* read (and parse) the first line (starting from p) that contains at
       least one character */
    if (!parse_oneline_sequence(fd, msa, p, seqno, 0, &aln_len, &error)) break;

    ++seqno;

    if (seqno == msa->count) break;
  }

  /* was the last block of sequences non-aligned? */
  if (error)
  {
    corax_msa_destroy(msa);
    return NULL;
  }

  if (seqno != msa->count)
  {
    corax_set_error(CORAX_ERROR_PHYLIP_SYNTAX,
                    "Found %d sequence(s) but expected %d",
                    seqno,
                    msa->count);
    corax_msa_destroy(msa);
    return NULL;
  }

  /* update the length of the alignment read so far, which will be used as the
     offset when appending data to the end of the sequences */
  sumlen = aln_len;

  /* now read the remaining blocks */
  seqno           = 0;
  aln_len         = 0;
  int block_count = 2;
  while (1)
  {
    char *p = getnextline(fd);

    /* read (and parse) the first line (starting from p) that contains at
       least one character */
    if (!parse_oneline_sequence(fd, msa, p, seqno, sumlen, &aln_len, &error))
      break;

    seqno = (seqno + 1) % msa->count;

    /* if data for all sequences were read, then append the alignment length
       to the sum, and go for the next block */
    if (!seqno)
    {
      sumlen += aln_len;
      aln_len = 0;
      block_count++;
    }
  }

  /* was the last block of sequences non-aligned? */
  if (error)
  {
    corax_msa_destroy(msa);
    return NULL;
  }

  /* if seqno != 0 then there were more (or less) sequences than expected */
  if (seqno)
  {
    corax_set_error(CORAX_ERROR_PHYLIP_SYNTAX,
                    "Found %d sequences in block %d but expected %d",
                    seqno,
                    block_count,
                    msa->count);
    corax_msa_destroy(msa);
    return NULL;
  }
  if (sumlen != msa->length)
  {
    corax_set_error(CORAX_ERROR_PHYLIP_SYNTAX,
                    "Sequence length is %d but expected %d",
                    sumlen,
                    msa->length);
    corax_msa_destroy(msa);
    return NULL;
  }

  return msa;
}

CORAX_EXPORT corax_msa_t *corax_phylip_parse_sequential(corax_phylip_t *fd)
{
  int  i, j;
  long headerlen;

  corax_msa_t *msa = (corax_msa_t *)malloc(sizeof(corax_msa_t));
  if (!msa)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return NULL;
  }

  /* read header */
  if (!parse_header(
          fd->line, &(msa->count), &(msa->length), CORAX_PHYLIP_SEQUENTIAL))
  {
    free(msa);
    return NULL;
  }

  msa->sequence = (char **)calloc((size_t)(msa->count), sizeof(char *));
  msa->label    = (char **)calloc((size_t)(msa->count), sizeof(char *));
  if (!msa->label || !msa->sequence)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    corax_msa_destroy(msa);
    return NULL;
  }

  for (i = 0; i < msa->count; ++i)
  {
    msa->sequence[i] = (char *)malloc((size_t)(msa->length + 1) * sizeof(char));
    if (!msa->sequence[i])
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Unable to allocate enough memory.");
      corax_msa_destroy(msa);
      return NULL;
    }
    msa->sequence[i][msa->length] = 0;
  }

  /* assume the worst: whole sequence in one line + long label. *
   * Doing reallocate+copy really hurts on long alignments,     *
   * so we better waste some memory here to avoid reallocation. */
  reallocline(fd, (size_t)msa->length + 300);

  /* read sequences */
  int seqno = 0;
  while (1)
  {
    /* get next line */
    fd->line = getnextline(fd);
    char *p  = fd->line;

    /* if no more lines break */
    if (!p) break;

    /* skip whitespace before sequence header */
    while (*p && whitespace(*p)) ++p;

    /* restart loop if blank line */
    if (!*p) continue;

    /* error if there are more sequences than specified */
    if (seqno == msa->count)
    {
      corax_set_error(CORAX_ERROR_PHYLIP_SYNTAX,
                      "Found at least %d sequences but expected %d",
                      seqno + 1,
                      msa->count);
      corax_msa_destroy(msa);
      return NULL;
    }

    /* find first delimiter after header */
    headerlen = get_headerlen(p);

    /* headerlen cannot be zero */
    assert(headerlen > 0);

    /* store sequence header */
    msa->label[seqno] = (char *)malloc((size_t)(headerlen + 1) * sizeof(char));
    if (!msa->label[seqno])
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Unable to allocate enough memory.");
      corax_msa_destroy(msa);
      return NULL;
    }
    memcpy(msa->label[seqno], p, (size_t)headerlen);
    msa->label[seqno][headerlen] = 0;

    p += headerlen;

    /* go through possibly multiple sequence data lines */
    j = 0;
    while (1)
    {
      /* read sequence data */
      int chars_count = dfa_parse(fd, msa, p, seqno, j);
      if (chars_count == -1)
      {
        corax_msa_destroy(msa);
        return NULL;
      }

      j += chars_count;

      /* break if we read all sequence data */
      if (j == msa->length) break;

      p = getnextline(fd);

      if (!p)
      {
        corax_set_error(
            CORAX_ERROR_PHYLIP_SYNTAX,
            "Sequence %d (%.100s) has %d characters but expected %d",
            seqno + 1,
            msa->label[seqno],
            j,
            msa->length);
        corax_msa_destroy(msa);
        return NULL;
      }
    }

    ++seqno;
  }

  if (seqno != msa->count)
  {
    corax_set_error(CORAX_ERROR_PHYLIP_SYNTAX,
                    "Found %d sequence(s) but expected %d",
                    seqno,
                    msa->count);
    corax_msa_destroy(msa);
    return NULL;
  }

  return msa;
}

corax_msa_t *corax_phylip_load(const char *fname, corax_bool_t interleaved)
{
  corax_phylip_t *fd = corax_phylip_open(fname, corax_map_generic);
  if (!fd) return NULL;

  corax_msa_t *msa = interleaved ? corax_phylip_parse_interleaved(fd)
                                 : corax_phylip_parse_sequential(fd);

  corax_phylip_close(fd);

  return msa;
}

CORAX_EXPORT int corax_phylip_save(const char *       out_fname,
                                   const corax_msa_t *msa)
{
  if (!msa)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "MSA structure is NULL");
    return CORAX_FAILURE;
  }

  if (!out_fname)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "File name (out_fname) is NULL");
    return CORAX_FAILURE;
  }

  FILE *f = fopen(out_fname, "w");

  if (!f)
  {
    corax_set_error(CORAX_ERROR_FILE_OPEN, "Cannot open file: %s", out_fname);
    return CORAX_FAILURE;
  }

  fprintf(
      f, "%lu %lu\n", (unsigned long)msa->count, (unsigned long)msa->length);

  unsigned long i;
  for (i = 0; i < (unsigned long)msa->count; ++i)
  {
    fprintf(f, "%s    %s\n", msa->label[i], msa->sequence[i]);
  }

  fclose(f);

  return CORAX_SUCCESS;
}

CORAX_EXPORT void corax_msa_destroy(corax_msa_t *msa)
{
  if (!msa) return;

  int i;

  if (msa->label)
  {
    for (i = 0; i < msa->count; ++i)
      if (msa->label[i]) free(msa->label[i]);
    free(msa->label);
  }

  if (msa->sequence)
  {
    for (i = 0; i < msa->count; ++i)
      if (msa->sequence[i]) free(msa->sequence[i]);
    free(msa->sequence);
  }

  free(msa);
}
