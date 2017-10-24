from sys import path

import isolib as il

il.padova_interpolated_isomake(['iphas', '2MASS'],
                               {'gR': 'r_INT', 'gI': 'i_INT', 'Ha': 'Ha_INT',
                                'J': 'J_2MASS', 'H': 'H_2MASS',
                                'Ks': 'Ks_2MASS'},
                               "padova_iphas-2MASS.txt",
                               ['r_INT', 'i_INT', 'Ha_INT', 'J_2MASS',
                                'H_2MASS', 'Ks_2MASS'])
