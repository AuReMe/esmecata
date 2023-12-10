import pandas as pd

from esmecata.analysis import normalise_dataframe

def test_normalisation_default_3_digits():
    input_dataframe = pd.DataFrame([[5, 5, 1, 2, 3],
                                    [15, 5, 1, 3, 2],
                                    [2, 3, 3, 10, 10],
                                    [2, 2, 4, 3, 3]],
                                    columns=['1.1.2.1', '1.1.2.4', '2.1.1.1', '3.1.1.1', '3.1.1.2'],
                                    index=['O1a', 'O1b', 'O2', 'O3'])

    taxa_name = {'T1': ['O1a', 'O1b'], 'T2': ['O2'], 'T3': ['O3']}

    nb_digit = 3

    normalised_dataframe = normalise_dataframe(input_dataframe, taxa_name, nb_digit)

    expected_dataframe = pd.DataFrame([[0.7692307692307693, 0.1282051282051282, 0.10256410256410256],
                                    [0.2222222222222222, 0.3333333333333333, 0.4444444444444444],
                                    [0.2777777777777778, 0.5555555555555556, 0.16666666666666666]],
                                    columns=['T1', 'T2', 'T3'],
                                    index=['1.1.2', '2.1.1', '3.1.1'])

    pd.testing.assert_frame_equal(normalised_dataframe, expected_dataframe)


def test_normalisation_default_4_digits():
    input_dataframe = pd.DataFrame([[5, 5, 1, 2, 3],
                                    [15, 5, 1, 3, 2],
                                    [2, 3, 3, 10, 10],
                                    [2, 2, 4, 3, 3]],
                                    columns=['1.1.2.1', '1.1.2.4', '2.1.1.1', '3.1.1.1', '3.1.1.2'],
                                    index=['O1a', 'O1b', 'O2', 'O3'])

    taxa_name = {'T1': ['O1a', 'O1b'], 'T2': ['O2'], 'T3': ['O3']}

    nb_digit = 4

    normalised_dataframe = normalise_dataframe(input_dataframe, taxa_name, nb_digit)

    expected_dataframe = pd.DataFrame([[0.833333, 0.083333, 0.083333],
                                    [0.666667, 0.200000, 0.133333],
                                    [0.222222, 0.333333, 0.444444],
                                    [0.277778, 0.555556, 0.166667],
                                    [0.277778, 0.555556, 0.166667]],
                                    columns=['T1', 'T2', 'T3'],
                                    index=['1.1.2.1', '1.1.2.4', '2.1.1.1', '3.1.1.1', '3.1.1.2'])

    pd.testing.assert_frame_equal(normalised_dataframe, expected_dataframe)
      
      
      
      