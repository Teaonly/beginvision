{
    'target_defaults': {
       'default_configuration': 'Debug',
        'configurations': {
            'Debug': {
                'include_dirs': [
                    '..',
                    '../third_party',
                ],
                'cflags': [
                    '-g',
                ],
            },
        },
    },
    'targets': [
    {
        'target_name': 'io',
        'type': 'executable',
        'sources': [
            'io.cpp',
        ],
    },
    {
        'target_name': 'filter',
        'type': 'executable',
        'sources': [
            'filter.cpp',
        ],
    },
    {
        'target_name': 'resize',
        'type': 'executable',
        'sources': [
            'resize.cpp',
        ],
    },
    {
        'target_name': 'harris',
        'type': 'executable',
        'sources': [
            'harris.cpp',
        ],
    },
    {
        'target_name': 'blob',
        'type': 'executable',
        'sources': [
            'blob.cpp',
        ],
    },
    {
        'target_name': 'lbfgs',
        'type': 'executable',
        'sources': [
            'lbfgs.cpp',
        ],
    },
    {
        'target_name': 'lkp',
        'type': 'executable',
        'sources': [
            'lkpTrack.cpp',
        ],
    },
    {
        'target_name': 'sift',
        'type': 'executable',
        'sources': [
            'sift.cpp',
        ],
    },
    {
        'target_name': 'affine',
        'type': 'executable',
        'sources': [
            'affineMatch.cpp',
        ],
    },

    ],
}
