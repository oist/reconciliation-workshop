include(targetSourcesLocal)

target_sources_local(
    corax-bench PRIVATE PLLTreeInfo.cpp PLLUnrootedTree.cpp LibpllParsers.cpp Model.cpp bench.cpp
)
