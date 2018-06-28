import os
import unittest

import pandas as pd
from firecloud import api as fapi

import kco


class CommandTests(unittest.TestCase):
    def test_create_sample_sheet_dir(self):
        kco.sample_sheet.main(['-d', 'inputs', '-o', 'test_sample_sheet.txt'])
        df = pd.read_table('test_sample_sheet.txt', header=None, names=['sample', 'path'])
        self.assertTrue(df.shape[1], 2)
        self.assertTrue(df.shape[0], 5)
        self.assertTrue(len(df['sample'].unique()) == 3)
        os.remove('test_sample_sheet.txt')

    def test_create_sample_sheet_r1_r2(self):
        kco.sample_sheet.main(['-d', 'inputs', '-o', 'test_sample_sheet.txt', '-f', 'r1_r2'])
        df = pd.read_table('test_sample_sheet.txt', header=None, names=['sample', 'r1', 'r2'])
        self.assertTrue(df.shape[1], 3)
        self.assertTrue(df.shape[0], 5)
        self.assertTrue(len(df['sample'].unique()) == 3)
        os.remove('test_sample_sheet.txt')

    def test_create_sample_sheet_r1_r2_i1(self):
        kco.sample_sheet.main(['-d', 'inputs', '-o', 'test_sample_sheet.txt', '-f', 'r1_r2_i1'])
        df = pd.read_table('test_sample_sheet.txt', header=None, names=['sample', 'r1', 'r2', 'i1'])
        self.assertTrue(df.shape[1], 4)
        self.assertTrue(df.shape[0], 5)
        self.assertTrue(len(df['sample'].unique()) == 3)
        os.remove('test_sample_sheet.txt')

    def test_upload_sample_sheet(self):
        kco.sample_sheet.main(['-d', 'inputs', '-o', 'test_sample_sheet.txt'])
        kco.upload.main(['--dry-run', '-w', 'regev-development/test', 'test_sample_sheet.txt'])
        os.remove('test_sample_sheet.txt')

    def test_run(self):
        kco.fc_run.main(['-m', 'broadgdac/echo', '-i', 'test.json', '-w', 'regev-development/test'])

    def test_inputs(self):
        kco.fc_inputs.main(['-m', 'broadgdac/echo', '-o', 'test.json'])
        os.remove('test.json')


if __name__ == '__main__':
    unittest.main()
