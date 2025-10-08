include(FetchContent)

FetchContent_Declare(
            CLI11
            GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
            GIT_TAG v2.3.2
            SYSTEM
)

FetchContent_Declare(
      tomlplusplus
      GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
      GIT_TAG        v3.4.0 # or latest
)


FetchContent_Declare(
        kamping
        GIT_REPOSITORY https://github.com/kamping-site/kamping.git
        GIT_TAG e2c5e6e
        SYSTEM
)

FetchContent_Declare(
        malloc_count
        GIT_REPOSITORY https://github.com/bingmann/malloc_count.git
        GIT_TAG ddd0565
        SYSTEM
)

FetchContent_Declare(
        absl
        GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
        GIT_TAG 4447c75
        SYSTEM
)

FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 10.1.1
        SYSTEM
)

FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 10.1.1
        SYSTEM
)
