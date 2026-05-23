import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer


def preprocess(df: pd.DataFrame):
    """Preprocess the clean twins data into a model-ready design matrix.

    Drops non-covariate columns, paternal variables (severe missingness), and
    high-cardinality geographic identifiers. Imputes missing values (median for
    ordinal, most frequent for nominal), one-hot encodes nominal categoricals,
    and constructs an explicit pair identifier for cluster-based inference.

    Returns:
        X (pd.DataFrame): model-ready covariate matrix (imputed + encoded)
        T (pd.Series): treatment indicator (1 = heavier twin)
        Y (pd.Series): outcome (1 = first-year mortality)
        exposure (pd.Series): any_exposure indicator for subgroup analysis
        pair_id (pd.Series): pair identifier for cluster bootstrap
    """
    # Step 1: Reset index and build pair_id
    df = df.reset_index(drop=True)
    pair_id = df.index // 2

    # Step 2: Pull out treatment, outcome, exposure
    T = df['treatment']
    Y = df['outcome']
    exposure = df['any_exposure']

    # Step 3: Define column groups
    ordinal_cols = ['mager8', 'meduc6', 'mpre5', 'adequacy', 'gestat10',
                    'cigar6', 'drink5', 'nprevistq',
                    'dlivord_min', 'dtotord_min']

    nominal_cols = ['pldel', 'birattnd', 'ormoth', 'mrace', 'dmar',
                    'birmon', 'csex', 'anemia', 'cardiac', 'lung', 'diabetes',
                    'herpes', 'hydra', 'hemo', 'chyper', 'phyper', 'eclamp',
                    'incervix', 'pre4000', 'preterm', 'renal', 'rh', 'uterine',
                    'othermr', 'tobacco', 'alcohol', 'crace', 'data_year',
                    'bord', 'brstate_reg', 'stoccfipb_reg', 'mplbir_reg']

    # Step 4: Impute ordinal columns with median
    ord_imputer = SimpleImputer(strategy='median')
    ordinal_imputed = pd.DataFrame(ord_imputer.fit_transform(df[ordinal_cols]), columns=ordinal_cols)

    # Step 5: Impute nominal columns with mode
    nom_imputer = SimpleImputer(strategy='most_frequent')
    nominal_imputed = pd.DataFrame(nom_imputer.fit_transform(df[nominal_cols]), columns=nominal_cols)

    # Step 6: One-hot encode the nominal columns
    nominal_encoded = pd.get_dummies(nominal_imputed, columns=nominal_cols,
        drop_first=True
    )

    # Step 7: Combine ordinal (numeric) + encoded nominal into final X
    X = pd.concat([ordinal_imputed, nominal_encoded], axis=1)
    X = X.astype(float)
    
    return X, T, Y, exposure, pair_id