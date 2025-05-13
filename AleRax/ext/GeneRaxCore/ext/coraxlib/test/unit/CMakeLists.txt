set(CORAX_TEST_DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data)

# Register the unittests 
corax_register_test(test_dummy FILES src/dummy.cpp)
corax_register_test(test_corax_assert FILES src/corax_assert.cpp)
corax_register_test(test_corax_assert_exceptions FILES src/corax_assert.cpp DEFINES CORAX_ASSERT_USE_EXCEPTIONS)
corax_register_test(test_newick FILES src/newick.cpp DEFINES DATAPATH=${CORAX_TEST_DATA_DIR})
corax_register_test(test_msa FILES src/msa.cpp DEFINES DATAPATH=${CORAX_TEST_DATA_DIR})
corax_register_test(test_evolmodel FILES src/evolmodel.cpp)
corax_register_test(test_nni FILES src/nni.cpp)
corax_register_test(test_partition FILES src/partition.cpp)
corax_register_test(test_absdiff FILES src/absdiff.cpp)
