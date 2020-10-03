# Copyright (c) 2020 Google LLC
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""Builds a table with SMILES and yield information."""

import collections
import dataclasses

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def location_to_row_col(location, block, plate):
    """Converts a block location to (row, col) on the plate.

    Args:
        location: Text location in the block, e.g. "1:A".
        block: Integer block number.
        plate: Integer plate number.

    Returns:
        Tuple of (row, col) integers; the location on the plate.
    """
    if plate == 3:
        row_letter, col = location.split(':')
    else:
        col, row_letter = location.split(':')
    col = int(col)
    row = ord(row_letter) - 64
    if block == 2:
        col += 24
    elif block == 3:
        row += 16
    elif block == 4:
        row += 16
        col += 24
    return row, col


def read_yield_data():
    """Reads location/yield data from the yield_data/ directory.

    Returns:
        DataFrame with the following columns:
            * plate
            * row
            * col
            * yield
    """
    data = []
    for plate in [1, 2, 3]:
        for block in [1, 2, 3, 4]:
            filename = f'yield_data/plate{plate}.{block}.csv'
            print(filename)
            df = pd.read_csv(filename)
            mask = df.Location.isna()
            if mask.any():
                print(df[mask])
                df = df[~mask]
            locations = df.apply(
                lambda x: location_to_row_col(x.Location, block, plate), axis=1,
                result_type='expand')
            locations.rename(columns={0: 'row', 1: 'col'}, inplace=True)
            df = pd.concat([df, locations], axis=1)
            df['plate'] = plate
            df = df[['plate', 'row', 'col', 'product_scaled']]
            df.rename(columns={'product_scaled': 'yield'}, inplace=True)
            data.append(df)
    return pd.concat(data, ignore_index=True)


def read_compound_data():
    """Reads location/compound data from the layout/ and smiles/ directories.

    Returns:
        DataFrame with the following columns:
            * plate
            * row
            * col
            * additive
            * additive_number
            * additive_smiles
            * aryl_halide
            * aryl_halide_number
            * aryl_halide_smiles
            * base
            * base_cas_number
            * base_smiles
            * ligand
            * ligand_cas_number
            * ligand_smiles
            * product_smiles
    """
    rows = pd.read_csv('layout/Table_S1.csv')
    cols = pd.read_csv('layout/Table_S2.csv')
    compounds = read_compounds()
    data = []
    for plate in [1, 2, 3]:
        for row in range(1, 33):
            row_data = rows[rows.Row == row].iloc[0]
            for col in range(1, 49):
                col_data = cols[cols.Column == col].iloc[0]
                base = compounds['base'][col_data.Base]
                ligand = compounds['ligand'][row_data.Ligand]
                record = {
                    'plate': plate,
                    'row': row,
                    'col': col,
                    'base': col_data.Base,
                    'base_cas_number': base.cas_number,
                    'base_smiles': base.smiles,
                    'ligand': row_data.Ligand,
                    'ligand_cas_number': ligand.cas_number,
                    'ligand_smiles': ligand.smiles,
                }
                if pd.notnull(row_data[f'Additive (Plate {plate})']):
                    record['additive_number'] = (
                        int(row_data[f'Additive (Plate {plate})']))
                    additive = compounds['additive'][record['additive_number']]
                    record['additive'] = additive.name
                    record['additive_smiles'] = additive.smiles
                if pd.notnull(col_data['Aryl Halide']):
                    record['aryl_halide_number'] = int(col_data['Aryl Halide'])
                    aryl_halide = (
                        compounds['aryl_halide'][record['aryl_halide_number']])
                    record['aryl_halide'] = aryl_halide.name
                    record['aryl_halide_smiles'] = aryl_halide.smiles
                data.append(record)
    df = pd.DataFrame(data)
    df['product_smiles'] = df.aryl_halide_smiles.map(get_product_smiles,
                                                     na_action='ignore')
    return df


@dataclasses.dataclass
class CompoundInfo:
    """Container for compound identifiers."""
    name: str
    smiles: str = ''
    cas_number: str = ''


def read_compounds():
    """Reads compound names and SMILES from the smiles/ directory.

    Returns:
        Dict of dicts mapping column -> number/name -> CompoundInfo.
    """
    data = collections.defaultdict(dict)
    for row in pd.read_csv('smiles/additive-list.csv').itertuples():
        data['additive'][row.component] = CompoundInfo(
            name=row.name, smiles=row.Additive_SMILES)
    for row in pd.read_csv('smiles/aryl_halide-list.csv').itertuples():
        data['aryl_halide'][row.component] = CompoundInfo(
            name=row.name, smiles=row.Aryl_halide_SMILES)
    for row in pd.read_csv('smiles/base-list.csv').itertuples():
        data['base'][row.name] = CompoundInfo(
            name=row.name, smiles=row.Base_SMILES, cas_number=row.CAS)
    for row in pd.read_csv('smiles/ligand-list.csv').itertuples():
        data['ligand'][row.name] = CompoundInfo(
            name=row.name, smiles=row.Ligand_SMILES, cas_number=row.CAS)
    return data


def get_product_smiles(aryl_halide):
    """Adds product SMILES to the DataFrame.

    Copied from https://github.com/Open-Reaction-Database/ord-schema/pull/88.

    Args:
        aryl_halide: Text aryl halide SMILES.

    Returns:
         Product SMILES.
    """
    replacement = Chem.MolFromSmiles('NC1=CC=C(C)C=C1')
    query = Chem.MolFromSmarts('[Cl,Br,I]')
    molecule = Chem.MolFromSmiles(aryl_halide)
    products = AllChem.ReplaceSubstructs(molecule, query, replacement)
    assert len(products) == 1
    return Chem.MolToSmiles(products[0])


def main():
    yield_data = read_yield_data()
    compound_data = read_compound_data()
    df = compound_data.merge(yield_data,
                             on=['plate', 'row', 'col'],
                             validate='one_to_one')
    print(df.shape)
    df.to_csv('data_table.csv', index=False)


if __name__ == '__main__':
    main()
