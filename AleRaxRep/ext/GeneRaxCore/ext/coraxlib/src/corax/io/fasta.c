/*
    Copyright (C) 2015 Tomas Flouri

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

#define MEMCHUNK 4096

/* please note that these functions will return a pointer to a buffer
   allocated here for the query header and sequence. This buffers will
   be overwritten on the next call of query_getnext. */

/* define strchrnul in case this is not a GNU system */
static char *xstrchrnul(char *s, int c)
{
  char *r = strchr(s, c);
  if (r) return r;

  return (char *)s + strlen(s);
}

CORAX_EXPORT corax_fasta_t *corax_fasta_open(const char *        filename,
                                             const unsigned int *map)
{
  int            i;
  corax_fasta_t *fd = (corax_fasta_t *)malloc(sizeof(corax_fasta_t));
  if (!fd)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return NULL;
  }

  /* allocate space */

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

  fd->line[0] = 0;
  if (!fgets(fd->line, CORAX_LINEALLOC, fd->fp))
  {
    corax_set_error(
        CORAX_ERROR_FILE_SEEK, "Unable to read file (%s)", filename);
    fclose(fd->fp);
    free(fd);
    return NULL;
  }
  fd->lineno = 1;

  return fd;
}

CORAX_EXPORT int corax_fasta_rewind(corax_fasta_t *fd)
{
  int i;

  rewind(fd->fp);

  /* reset stripped char frequencies */
  fd->stripped_count = 0;
  for (i = 0; i < 256; i++) fd->stripped[i] = 0;

  fd->line[0] = 0;
  if (!fgets(fd->line, CORAX_LINEALLOC, fd->fp))
  {
    corax_set_error(CORAX_ERROR_FILE_SEEK, "Unable to rewind and cache data");
    return CORAX_FAILURE;
  }
  fd->lineno = 1;

  return CORAX_SUCCESS;
}

CORAX_EXPORT void corax_fasta_close(corax_fasta_t *fd)
{
  fclose(fd->fp);
  free(fd);
}

CORAX_EXPORT int corax_fasta_getnext(corax_fasta_t *fd,
                                     char **        head,
                                     long *         head_len,
                                     char **        seq,
                                     long *         seq_len,
                                     long *         seqno)
{
  void *mem;
  long  head_alloc = MEMCHUNK;
  long  seq_alloc  = MEMCHUNK;

  *head_len = 0;
  *seq_len  = 0;

  /* allocate sequence buffers */
  *head = (char *)malloc((size_t)(head_alloc));
  if (!(*head))
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return CORAX_FAILURE;
  }

  *seq = (char *)malloc((size_t)(seq_alloc));
  if (!(*seq))
  {
    free(*head);
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return CORAX_FAILURE;
  }

  /* read line and increase line number */

  while (fd->line[0])
  {
    /* read header */

    if (fd->line[0] != '>')
    {
      corax_set_error(CORAX_ERROR_FASTA_INVALIDHEADER,
                      "Illegal header line in query fasta file");
      free(*head);
      free(*seq);
      return CORAX_FAILURE;
    }

    long headerlen;
    if (strchr(fd->line + 1, '\r'))
      headerlen = xstrchrnul(fd->line + 1, '\r') - (fd->line + 1);
    else
      headerlen = xstrchrnul(fd->line + 1, '\n') - (fd->line + 1);

    *head_len = headerlen;

    if (headerlen + 1 > head_alloc)
    {
      head_alloc = headerlen + 1;
      mem        = realloc(*head, (size_t)(head_alloc));
      if (!mem)
      {
        corax_set_error(CORAX_ERROR_MEM_ALLOC,
                        "Unable to allocate enough memory.");
        free(*head);
        free(*seq);
        return CORAX_FAILURE;
      }
      *head = (char *)mem;
    }

    memcpy(*head, fd->line + 1, (size_t)headerlen);
    *(*head + headerlen) = 0;

    /* get next line */

    fd->line[0] = 0;
    if (!fgets(fd->line, CORAX_LINEALLOC, fd->fp))
    { /* do nothing */
    }
    fd->lineno++;

    /* read sequence */

    *seq_len = 0;

    while (fd->line[0] && (fd->line[0] != '>'))
    {
      char  c;
      char  m;
      char *p = fd->line;

      while ((c = *p++))
      {
        m = (char)fd->chrstatus[(int)c];
        switch (m)
        {
        case 0:
          /* character to be stripped */
          fd->stripped_count++;
          fd->stripped[(int)c]++;
          break;

        case 1:
          /* legal character */
          if (*seq_len + 1 > seq_alloc)
          {
            seq_alloc += MEMCHUNK;
            mem = realloc(*seq, (size_t)(seq_alloc));
            if (!mem)
            {
              corax_set_error(CORAX_ERROR_MEM_ALLOC,
                              "Unable to allocate enough memory.");
              free(*head);
              free(*seq);
              return CORAX_FAILURE;
            }
            *seq = (char *)mem;
          }
          *(*seq + *seq_len) = c;
          (*seq_len)++;

          break;

        case 2:
          /* fatal character */
          if (c >= 32)
          {
            corax_set_error(CORAX_ERROR_FASTA_ILLEGALCHAR,
                            "illegal character '%c' "
                            "on line %ld in the fasta file",
                            c,
                            fd->lineno);
          }
          else
          {
            corax_set_error(CORAX_ERROR_FASTA_UNPRINTABLECHAR,
                            "illegal unprintable character "
                            "%#.2x (hexadecimal) on line %ld "
                            "in the fasta file",
                            c,
                            fd->lineno);
          }
          return CORAX_FAILURE;

        case 3:
          /* silently stripped chars */
          break;
        }
      }

      fd->line[0] = 0;
      if (!fgets(fd->line, CORAX_LINEALLOC, fd->fp))
      { /* do nothing */
      }
      fd->lineno++;
    }

    /* add zero after sequence */

    if (*seq_len + 1 > seq_alloc)
    {
      seq_alloc += MEMCHUNK;
      mem = realloc(*seq, (size_t)seq_alloc);
      if (!mem)
      {
        corax_set_error(CORAX_ERROR_MEM_ALLOC,
                        "Unable to allocate enough memory.");
        free(*head);
        free(*seq);
        return CORAX_FAILURE;
      }
      *seq = (char *)mem;
    }
    *(*seq + *seq_len) = 0;

    fd->no++;
    *seqno = fd->no;

    return CORAX_SUCCESS;
  }

  corax_set_error(CORAX_ERROR_FILE_EOF, "End of file\n");
  free(*head);
  free(*seq);
  return CORAX_FAILURE;
}

CORAX_EXPORT long corax_fasta_getfilesize(const corax_fasta_t *fd)
{
  return fd->filesize;
}

CORAX_EXPORT long corax_fasta_getfilepos(const corax_fasta_t *fd)
{
  return ftell(fd->fp);
}

CORAX_EXPORT corax_msa_t *corax_fasta_load(const char *fname)
{
  int i;

  char *seq = NULL;
  char *hdr = NULL;
  long  seqlen;
  long  hdrlen;
  long  seqno;

  long   alloc_count = 0;
  long   alloc_chunk = 0;
  size_t alloc_size  = 0;

  corax_fasta_t *fp = corax_fasta_open(fname, corax_map_generic);
  if (!fp)
  {
    assert(corax_errno);
    return NULL;
  }

  corax_msa_t *msa = (corax_msa_t *)calloc(1, sizeof(corax_msa_t));
  if (!msa)
  {
    corax_fasta_close(fp);
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return NULL;
  }

  /* read FASTA sequences and make sure they are all of the same length */
  msa->length = -1;
  msa->count  = 0;
  for (i = 0; corax_fasta_getnext(fp, &hdr, &hdrlen, &seq, &seqlen, &seqno);
       ++i)
  {
    if (msa->length == -1)
    {
      msa->length = (int)seqlen;

      /* estimate number of sequences from file size and sequence length
       * (it's better to overestimate slightly to avoid reallocations) */
      alloc_chunk = (long)ceil(1.1 * fp->filesize / (seqlen + hdrlen));
    }
    else if (msa->length != seqlen)
    {
      msa->count = i;
      corax_msa_destroy(msa);
      corax_fasta_close(fp);
      corax_set_error(CORAX_ERROR_FASTA_NONALIGNED,
                      "FASTA file does not contain equal size sequences: "
                      "sequence %d has length of %ld (expected: %d)",
                      i,
                      seqlen,
                      msa->length);
      return NULL;
    }

    if (i >= alloc_count)
    {
      /* allocate arrays to store FASTA headers and sequences */
      alloc_count += alloc_chunk;
      alloc_size    = (size_t)alloc_count * sizeof(char *);
      msa->label    = (char **)realloc(msa->label, alloc_size);
      msa->sequence = (char **)realloc(msa->sequence, alloc_size);
    }

    msa->label[i]    = hdr;
    msa->sequence[i] = seq;
  }

  msa->count = i;

  /* trim label and sequence arrays to their actual size */
  alloc_size    = (size_t)msa->count * sizeof(char *);
  msa->label    = (char **)realloc(msa->label, alloc_size);
  msa->sequence = (char **)realloc(msa->sequence, alloc_size);

  /* close FASTA file */
  corax_fasta_close(fp);

  /* did we stop reading the file because we reached EOF? */
  if (corax_errno != CORAX_ERROR_FILE_EOF)
  {
    corax_msa_destroy(msa);
    return NULL;
  }

  corax_errno = 0;

  return msa;
}
