pipeline {
    agent {
        label 'cme-eastwatch'
    }
    stages {
        stage('BuildAndTest') {
            matrix {
                axes {
                    axis {
                        name 'DOCKERFILE'
                        values 'dockerfile-gcc-7', 'dockerfile-gcc-8', 'dockerfile-gcc-9', 'dockerfile-gcc-10', 'dockerfile-gcc-11', 'dockerfile-gcc-12', 'dockerfile-clang'
                    }
                    axis {
                        name 'RELEASE_TYPE'
                        values 'debug', 'release'
                    }
                    axis {
                        name 'BUILD_DIFFICULTY'
                        values '0', '1'
                    }
                }
                excludes {
                    exclude {
                        axis {
                            name 'DOCKERFILE'
                            values 'dockerfile-gcc-7', 'dockerfile-gcc-8', 'dockerfile-gcc-9', 'dockerfile-gcc-10'
                        }
                        axis {
                            name 'RELEASE_TYPE'
                            values 'release'
                        }
                        axis {
                            name 'BUILD_DIFFICULTY'
                            values '0'
                        }
                    }
                }
                agent {
                    dockerfile {
                        reuseNode true
                        filename "${DOCKERFILE}"
                        dir 'ci'
                    }
                }
                environment {
                    BUILD_DIR = "BUILD_DIR_${DOCKERFILE}_${RELEASE_TYPE}_${BUILD_DIFFICULTY}"
                }
                stages {
                    stage('Submodules') {
                        steps {
                            catchError(buildResult: 'FAILURE', stageResult: 'FAILURE') {
                                sh """#!/bin/bash
                                echo Get submodules for ${DOCKERFILE} - ${RELEASE_TYPE} - DIFFICULTY_PRED=${BUILD_DIFFICULTY}
                                set -exo pipefail
                                git submodule update --recursive --init
                                """
                            }
                        }
                    }
                    stage('Build') {
                        steps {
                            catchError(buildResult: 'FAILURE', stageResult: 'FAILURE') {
                                sh """#!/bin/bash
                                echo Do Build for ${DOCKERFILE} - ${RELEASE_TYPE} - DIFFICULTY_PRED=${BUILD_DIFFICULTY}
                                set -exo pipefail
                                rm -fr ${BUILD_DIR} && mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
                                cmake -DCMAKE_BUILD_TYPE=${RELEASE_TYPE} -DCORAX_BUILD_TESTS=1 -DCORAX_BUILD_DIFFICULTY_PREDICTION=${BUILD_DIFFICULTY} .. 2>&1 |tee cmake.out
                                make 2>&1 |tee make.out"""
                            }
                        }
                    }
                    stage('Test') {
                        steps {
                            catchError(buildResult: 'FAILURE', stageResult: 'FAILURE') {
                                sh """#!/bin/bash
                                echo Run Tests for ${DOCKERFILE} - ${RELEASE_TYPE} - DIFFICULTY_PRED=${BUILD_DIFFICULTY}
                                set -exo pipefail
                                cd ${BUILD_DIR}/test && make test"""
                            }
                        }
                    }
                }
            }
        }
        stage('Benchmark') {
            agent {
                dockerfile {
                    reuseNode true
                    filename 'dockerfile-gcc-11-bench'
                    dir 'ci'
                }
            }
            environment {
                BUILD_DIR = "BUILD_DIR_benchmark" // for now, only use gcc-11 to build and run the benchmark
            }
            stages {
                stage('Build Benchmark') {
                    steps {
                        catchError(buildResult: 'FAILURE', stageResult: 'FAILURE') {
                            sh """#!/bin/bash
                            echo Do Build for Benchmark
                            set -exo pipefail
                            rm -fr ${BUILD_DIR} && mkdir -p ${BUILD_DIR} && cd ${BUILD_DIR}
                            cmake -DCORAX_BUILD_BENCHMARKS=1 .. 2>&1 |tee cmake.out
                            make 2>&1 |tee make.out"""
                        }
                    }
                }
                stage('Run Benchmark') {
                    steps {
                        catchError(buildResult: 'FAILURE', stageResult: 'FAILURE') {
                            sh """#!/bin/bash
                            echo Run Benchmarks
                            set -exo pipefail
                            cd ${BUILD_DIR}
                            ./bench/corax-bench 2>&1 |tee bench.out"""
                        }
                    }
                }
            }
        }
        stage('Softwipe') {
            agent {
                dockerfile {
                    reuseNode true
                    filename 'dockerfile-softwipe'
                    dir 'ci'
                }
            }
            stages {
                stage('Run Softwipe') {
                    steps {
                        catchError(buildResult: 'FAILURE', stageResult: 'FAILURE') {
                            sh 'softwipe.py -CM . 2>&1 |tee softwipe_general.txt'
                        }
                    }
                }
            }
        }
    }
}