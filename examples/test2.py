from sys import path

import isolib as il

il.padova_interpolated_isomake(['iphas2', 'UKIDSS'],
                               {'gR': 'r_INT', 'gI': 'i_INT', 'Ha': 'Ha_INT',
                                'J': 'J_UKIDSS', 'H': 'H_UKIDSS',
                                'K': 'K_UKIDSS'},
                               "padova_iphas-UKIDSS.txt",
                               ['r_INT', 'i_INT', 'Ha_INT', 'J_UKIDSS',
                                'H_UKIDSS', 'K_UKIDSS'])
