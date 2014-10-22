from sldb.common.mutations import *

def _same_dict(d1, d2):
    return cmp(d1, d2) == 0

def test_germline_comparison():
    germline = '''
    CAGGTTCAGCTGGTGCAGTCTGGAGCT---GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCT
    TCTGGTTACACCTTT------------ACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTT
    GAGTGGATGGGATGGATCAGCGCTTAC------AATGGTAACACAAACTATGCACAGAAGCTCCAG---GGCAGA
    GTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCC
    GTGTATTACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGCCAGGGCACCCTGGTCACCGTC
    TCCTCAG
    '''.replace('\n', '').strip()

    seq = '''
    NNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGCAAGGCT
    TCTGGAGGCACCTTC------------AGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTT
    GAGTGGATGGGATGGATCAGCGCTTAC------AATGGTAACACAAACTATGCACAGAAGCTCCAG---GGCAGA
    GTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCC
    GTGTATTACTGTGCGTGAGCAGCAGCTGGTACGAATGGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTC
    TCCTCAG
    '''.replace('\n', '').strip()

    cdr3_len = 42

    mut = Mutations(germline, cdr3_len)
    mut.add_sequence(seq)

    regions, poss = mut.get_aggregate()

    assert _same_dict(regions,
        {  
            'all': {  
                'counts': {  
                    'total': {  
                        'sum':5,
                        'synonymous':0,
                        'nonconservative':4,
                        'conservative':1
                    },
                    'unique': {  
                        'sum':5,
                        'synonymous':0,
                        'nonconservative':4,
                        'conservative':1
                    }
                },
                'mutations': {  
                    'synonymous':[  

                    ],
                    'nonconservative':[  
                        (1,
                        (84,
                        'T',
                        'A'                )),
                        (1,
                        (85,
                        'T',
                        'G'                )),
                        (1,
                        (86,
                        'A',
                        'G'                )),
                        (1,
                        (375,
                        'C',
                        'A'                ))
                    ],
                    'conservative':[  
                        (1,
                        (116,
                        'G',
                        'C'                ))
                    ]
                }
            },
            'CDR3': {  
                'counts': {  
                    'total': {  
                        'sum':0,
                        'synonymous':0,
                        'nonconservative':0,
                        'conservative':0
                    },
                    'unique': {  
                        'sum':0,
                        'synonymous':0,
                        'nonconservative':0,
                        'conservative':0
                    }
                },
                'mutations': {  
                    'synonymous':[  

                    ],
                    'nonconservative':[  

                    ],
                    'conservative':[  

                    ]
                }
            },
            'CDR2': {  
                'counts': {  
                    'total': {  
                        'sum':0,
                        'synonymous':0,
                        'nonconservative':0,
                        'conservative':0
                    },
                    'unique': {  
                        'sum':0,
                        'synonymous':0,
                        'nonconservative':0,
                        'conservative':0
                    }
                },
                'mutations': {  
                    'synonymous':[  

                    ],
                    'nonconservative':[  

                    ],
                    'conservative':[  

                    ]
                }
            },
            'CDR1': {  
                'counts': {  
                    'total': {  
                        'sum':3,
                        'synonymous':0,
                        'nonconservative':3,
                        'conservative':0
                    },
                    'unique': {  
                        'sum':3,
                        'synonymous':0,
                        'nonconservative':3,
                        'conservative':0
                    }
                },
                'mutations': {  
                    'synonymous':[  

                    ],
                    'nonconservative':[  
                        (1,
                        (84,
                        'T',
                        'A'                )),
                        (1,
                        (85,
                        'T',
                        'G'                )),
                        (1,
                        (86,
                        'A',
                        'G'                ))
                    ],
                    'conservative':[  

                    ]
                }
            },
            'FR3': {  
                'counts': {  
                    'total': {  
                        'sum':0,
                        'synonymous':0,
                        'nonconservative':0,
                        'conservative':0
                    },
                    'unique': {  
                        'sum':0,
                        'synonymous':0,
                        'nonconservative':0,
                        'conservative':0
                    }
                },
                'mutations': {  
                    'synonymous':[  

                    ],
                    'nonconservative':[  

                    ],
                    'conservative':[  

                    ]
                }
            },
            'FR2': {  
                'counts': {  
                    'total': {  
                        'sum':1,
                        'synonymous':0,
                        'nonconservative':0,
                        'conservative':1
                    },
                    'unique': {  
                        'sum':1,
                        'synonymous':0,
                        'nonconservative':0,
                        'conservative':1
                    }
                },
                'mutations': {  
                    'synonymous':[  

                    ],
                    'nonconservative':[  

                    ],
                    'conservative':[  
                        (1,
                        (116,
                        'G',
                        'C'                ))
                    ]
                }
            },
            'FR1': {  
                'counts': {  
                    'total': {  
                        'sum':0,
                        'synonymous':0,
                        'nonconservative':0,
                        'conservative':0
                    },
                    'unique': {  
                        'sum':0,
                        'synonymous':0,
                        'nonconservative':0,
                        'conservative':0
                    }
                },
                'mutations': {  
                    'synonymous':[  

                    ],
                    'nonconservative':[  

                    ],
                    'conservative':[  

                    ]
                }
            },
            'FR4': {  
                'counts': {  
                    'total': {  
                        'sum':1,
                        'synonymous':0,
                        'nonconservative':1,
                        'conservative':0
                    },
                    'unique': {  
                        'sum':1,
                        'synonymous':0,
                        'nonconservative':1,
                        'conservative':0
                    }
                },
                'mutations': {  
                    'synonymous':[  

                    ],
                    'nonconservative':[  
                        (1,
                        (375,
                        'C',
                        'A'                ))
                    ],
                    'conservative':[  

                    ]
                }
            }
        })


    assert _same_dict(poss, 
        {  
            116: {  
                'synonymous':0,
                'nonconservative':0,
                'conservative':1
            },
            84: {  
                'synonymous':0,
                'nonconservative':1,
                'conservative':0
            },
            85: {  
                'synonymous':0,
                'nonconservative':1,
                'conservative':0
            },
            86: {  
                'synonymous':0,
                'nonconservative':1,
                'conservative':0
            },
            375: {  
                'synonymous':0,
                'nonconservative':1,
                'conservative':0
            }
        })
