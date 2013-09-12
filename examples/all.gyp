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
        'target_name': 'ex1',
        'type': 'executable',
        'sources': [
            'example1.cpp',
        ],
    },
    {
        'target_name': 'ex2',
        'type': 'executable',
        'sources': [
            'example2.cpp',
        ],
    },
    {
        'target_name': 'ex3',
        'type': 'executable',
        'sources': [
            'example3.cpp',
        ],
    },
    {
        'target_name': 'ex4',
        'type': 'executable',
        'sources': [
            'example4.cpp',
        ],
    },
    {
        'target_name': 'ex5',
        'type': 'executable',
        'sources': [
            'example5.cpp',
        ],
    },
    {
        'target_name': 'ex6',
        'type': 'executable',
        'sources': [
            'example6.cpp',
        ],
    },
    {
        'target_name': 'ex7',
        'type': 'executable',
        'sources': [
            'example7.cpp',
        ],
    },
    {
        'target_name': 'ex8',
        'type': 'executable',
        'sources': [
            'example8.cpp',
        ],
    },


    ],
}
