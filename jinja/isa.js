{
    "avx": {
        "fptypes": {
            "double": {"veclen": 4, "soalens": [2, 4]},
            "float": {"veclen": 8, "soalens": [4, 8]}
        },
        "extra_includes": ["immintrin.h", "qphix/avx/avx_utils.h"]
    },
    "avx2": {
        "fptypes": {
            "double": {"veclen": 4, "soalens": [2, 4]},
            "float": {"veclen": 8, "soalens": [4, 8]},
            "half": {"veclen": 8, "soalens": [4, 8]}
        },
        "extra_includes": ["immintrin.h", "qphix/avx/avx_utils.h"]
    },
    "avx512": {
        "fptypes": {
            "double": {"veclen": 8, "soalens": [4, 8]},
            "float": {"veclen": 16, "soalens": [4, 8, 16]},
            "half": {"veclen": 16, "soalens": [4, 8, 16]}
        },
        "extra_includes": ["immintrin.h"]
    },
    "mic": {
        "fptypes": {
            "double": {"veclen": 8, "soalens": [4, 8]},
            "float": {"veclen": 16, "soalens": [4, 8, 16]},
            "half": {"veclen": 16, "soalens": [4, 8, 16]}
        },
        "extra_includes": []
    },
    "qpx": {
        "fptypes": {
            "double": {"veclen": 4, "soalens": [4]}
        },
        "extra_includes": []
    },
    "scalar": {
        "fptypes": {
            "double": {"veclen": 1, "soalens": [1]},
            "float": {"veclen": 1, "soalens": [1]}
        },
        "extra_includes": []
    },
    "sse": {
        "fptypes": {
            "double": {"veclen": 2, "soalens": [2]},
            "float": {"veclen": 4, "soalens": [4]}
        },
        "extra_includes": ["immintrin.h", "qphix/sse/sse_utils.h"]
    }
}
