import pandas as pd
import numpy as np


def create_kz_gazetas_v_lte_0p4():
    t4 = '    '
    fname = 'kz_gazetas_v_lte_0p4'
    df = pd.read_csv(f'{fname}.csv')
    names = ['1p2', '4', '6', '10', '1000']
    a0s = np.linspace(0, 2, 10)
    fstr = f'def get_{fname}(a0, lob):\n'
    fstr += t4 + 'lobs = np.array([%s])\n' % ', '.join([xstr.replace('p', '.') for xstr in names])
    fstr += t4 + 'a0s = np.linspace(0, 2, 10)\n'
    fstr += t4 + 'kzs = np.array([\n'
    for name in names:

        kzs = np.interp(a0s, df[f'a0_lob_{name}'], df[f'kz_lob_{name}'])
        pstr = '[%s],' % ', '.join([f'{x:.4g}' for x in kzs])
        fstr += t4 + t4 + pstr + '\n'
    fstr += t4 + t4 + '])'

    print(fstr)


def create_kz_gazetas_v_gt_0p4():
    t4 = '    '
    fname = 'kz_gazetas_v_gt_0p4'
    df = pd.read_csv(f'{fname}.csv')
    names = ['1', '6', '1000']
    a0s = np.linspace(0, 2, 10)
    fstr = f'def get_{fname}(a0, lob):\n'
    fstr += t4 + 'lobs = np.array([%s])\n' % ', '.join([xstr.replace('p', '.') for xstr in names])
    fstr += t4 + 'a0s = np.linspace(0, 2, 10)\n'
    fstr += t4 + 'kzs = np.array([\n'
    for name in names:

        kzs = np.interp(a0s, df[f'a0_lob_{name}'], df[f'kz_lob_{name}'])
        pstr = '[%s],' % ', '.join([f'{x:.4g}' for x in kzs])
        fstr += t4 + t4 + pstr + '\n'
    fstr += t4 + t4 + '])'

    print(fstr)


def create_ky_gazetas():
    t4 = '    '
    fname = 'ky_gazetas'
    df = pd.read_csv(f'{fname}.csv')
    names = ['1', '2', '4', '6', '10', '1000']
    a0s = np.linspace(0, 2, 10)
    fstr = f'def get_{fname}(a0, lob):\n'
    fstr += t4 + 'lobs = np.array([%s])\n' % ', '.join([xstr.replace('p', '.') for xstr in names])
    fstr += t4 + 'a0s = np.linspace(0, 2, 10)\n'
    fstr += t4 + 'kys = np.array([\n'
    for name in names:

        kzs = np.interp(a0s, df[f'a0_lob_{name}'], df[f'ky_lob_{name}'])
        pstr = '[%s],' % ', '.join([f'{x:.4g}' for x in kzs])
        fstr += t4 + t4 + pstr + '\n'
    fstr += t4 + t4 + '])'

    print(fstr)


def gen_builder(fname, names, prefix):
    t4 = '    '

    df = pd.read_csv(f'{fname}.csv')

    a0s = np.linspace(0, 2, 10)
    fstr = f'def get_{fname}(a0, lob):\n'
    fstr += t4 + 'lobs = np.array([%s])\n' % ', '.join([xstr.replace('p', '.') for xstr in names])
    fstr += t4 + 'a0s = np.linspace(0, 2, 10)\n'
    fstr += t4 + f'{prefix}s = np.array([\n'
    for name in names:

        kzs = np.interp(a0s, df[f'a0_lob_{name}'], df[f'{prefix}_lob_{name}'])
        pstr = '[%s],' % ', '.join([f'{x:.4g}' for x in kzs])
        fstr += t4 + t4 + pstr + '\n'
    fstr += t4 + t4 + '])\n'
    fstr += t4 + f'return interp_twoways(a0, lob, a0s, lobs, {prefix}s)'

    print(fstr)


def create_cz_gazetas():
    fname = 'cz_gazetas_v_lte_0p4'
    names = ['1', '2', '4', '6', '10', '1000']
    prefix = 'cz'
    gen_builder(fname, names, prefix)


def create_czf_gazetas():
    fname = 'czf_gazetas_v_gt_0p4'
    names = ['1', '2', '4', '6', '1000']
    prefix = 'czf'
    gen_builder(fname, names, prefix)


def create_cy_gazetas_e_0p3():
    fname = 'cy_gazetas_v_e_0p3'
    names = ['1', '2', '4', '6', '10', '1000']
    prefix = 'cy'
    gen_builder(fname, names, prefix)


def create_cy_gazetas_e_0p5():
    fname = 'cy_gazetas_v_e_0p5'
    names = ['1', '2', '4', '6', '10', '1000']
    prefix = 'cy'
    gen_builder(fname, names, prefix)


def create_crx_gazetas():
    fname = 'crx_gazetas'
    names = ['1', '2', '4', '6', '10']
    prefix = 'crx'
    gen_builder(fname, names, prefix)


def create_cry_gazetas():
    fname = 'cry_gazetas'
    names = ['1', '2', '4', '6', '10', '1000']
    prefix = 'cry'
    gen_builder(fname, names, prefix)


def create_ct_gazetas():
    fname = 'ct_gazetas'
    names = ['1', '2', '3', '4', '6', '1000']
    prefix = 'ct'
    gen_builder(fname, names, prefix)



if __name__ == '__main__':
    # create_kz_gazetas_v_gt_0p4()
    # create_kz_gazetas_v_gt_0p4()
    create_crx_gazetas()
    create_cry_gazetas()
    create_ct_gazetas()