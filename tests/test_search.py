#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test the search module
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import pytest
import logging
import pysam
# import pandas

LOG = logging.getLogger(__name__)

try:
    GENOME_PATH = os.path.join(os.path.dirname(__file__), 'data/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz')
    assert os.path.exists(GENOME_PATH)
except AssertionError:
    LOG.error('Human reference genome file not found. Download the human reference genome')
    raise


@pytest.mark.parametrize(("pattern", 'mismatches', 'exphits_path', 'result_path'), [
    (b'TGGATGTGAAATGAGTCAAG', 3, 'data/TGGATGTGAAATGAGTCAAG-results.sam', 'Pattern_Matches_For_TGGATGTGAAATGAGTCAAG.txt'),
    (b'GGGTGGGGGGAGTTTGCTCC', 3, 'data/vegfa-site1-results.sam', 'Pattern_Matches_For_GGGTGGGGGGAGTTTGCTCC.txt'),
])
def test_search(pattern, mismatches, exphits_path, result_path):
    # TODO
    expected_hits = set()
    simplified_expected_hits = set()
    with open(exphits_path, 'rb') as exphits:
        for hit in exphits.readlines():
            # TODO use pysam to parse the expected result records if needed.
            expected_hits.add(hit)   
    
    # - Extract chromosome / start position / sequence
    for read in expected_hits:
        if read.startswith('@'):
            continue
        fields = read.split()
        hit = ','.join([fields[2], fields[3]])
        simplified_expected_hits.add(hit)

    predicted_hits = set()
    with open(result_path, 'rb') as predhits:
        for hit in predhits.readlines():
            fields = hit.split(',')
            hts = ','.join([fields[0], fields[2]])
            predicted_hits.add(hts)

    """
        Test might fail because of the following reasons
        - profile searching did not run completely through all trunks of all chromosomes
        - answer SAM file include unexpected matches with mismatches as high as 16 out of 20
    """
    assert (len(simplified_expected_hits.difference(predicted_hits)) == 0), "Expected all predicted!" 
    assert (len(predicted_hits.difference(simplified_expected_hits)) == 0), "Predicted all expected!" 
