#!/usr/bin/env python

"""
Project_Name: main, File_name: stage2_search.py
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 7/05/15 , Time:10:05 PM

Perform stage 2 in parallel
"""

import filters.constraints.looplengthConstraint as llc
import filters.noe.allNoes   as Noe
import filters.pcs.pcsfilter as Pfilter
import filters.rdc.rdcfilter as Rfilter
import filters.rmsd.RefRmsd as ref
import filters.rmsd.qcp as qcp
import filters.sequence.sequence_similarity as Sfilter
import ranking.RDCStageRank as rank
import utility.io_util as io
import utility.smotif_util as sm
import utility.stage2_util as uts2


def getfromDB(previous_smotif, current_ss, direction, database_cutoff):
    """
    :param database_cutoff:
    :param previous_smotif:
    :param current_ss:
    :param direction:
    :return:
    """
    psmotif = []
    for entry in previous_smotif:
        if 'smotif_def' == entry[0]:
            psmotif = entry[-1]

    if direction == 'left':
        previous_ss = psmotif[0]
    else:
        previous_ss = psmotif[1]

    # current_ss, previous_ss
    if direction == 'left':  # double check this implementation
        smotif_def = sm.getSmotif(current_ss, previous_ss)
    else:
        smotif_def = sm.getSmotif(previous_ss, current_ss)

    return sm.readSmotifDatabase(smotif_def, database_cutoff), smotif_def


def orderSSE(previous_smotif, current_sse, direction):
    """
    :param previous_smotif:
    :param current_sse:
    :param direction:
    :return:
    """

    previous_seq = []
    for entry in previous_smotif:
        if 'smotif_def' == entry[0]:
            previous_seq = entry

    ordered_SSE = []
    for sse_array in previous_seq[1]:
        ordered_SSE.append(sse_array)

    if direction == 'left':
        ordered_SSE.insert(0, current_sse)
    else:
        ordered_SSE.append(current_sse)
    return ordered_SSE


def SmotifSearch(task):
    """
    Main()
    :param task:
    :return:
    """
    index_array = task[0]
    stage = task[1]


    exp_data = io.readPickle("exp_data.pickle")
    exp_data_types = exp_data.keys()  # ['ss_seq', 'pcs_data', 'aa_seq', 'contacts']

    psmotif = uts2.getPreviousSmotif(index_array[0])

    current_ss, direction = uts2.getSS2(index_array[1])
    csmotif_data, smotif_def = getfromDB(psmotif, current_ss, direction, exp_data['database_cutoff'])

    if not csmotif_data:
        # If the smotif library doesn't exist
        # Terminate further execution
        return True

    """
    always narrow down to previous sse and current sse and operate on them individually
    """
    sse_ordered = orderSSE(psmotif, current_ss, direction)
    dump_log = []
    ref_rmsd = 0.0
    no_clashes = False

    # ************************************************************************************************
    # Main
    # The 'for' loop below iterates over all of the Smotifs and applies various filters
    # This is the place to add new filters as you desire. For starters, look at Sequence filter.
    # ************************************************************************************************

    # print "no of entries: ", len(csmotif_data)

    for i in range(0, len(csmotif_data)):

        # ************************************************
        # Applying different filters for the Smotif assembly
        # ************************************************

        # Exclude natives if needed
        tpdbid = csmotif_data[i][0][0]
        pdbid = tpdbid[0:4]

        if 'natives' in exp_data_types:
            natives = exp_data['natives']
            if pdbid in natives:
                continue
            # Stop further execution, but, iterate.
            else:
                pass

        if 'homologs' in exp_data_types:
            homologs = exp_data['homologs']
            if pdbid not in homologs:
                # Stop further execution, but, iterate.
                continue
            else:
                pass

        # ************************************************
        # RMSD filter using QCP method
        # quickly filters non-overlapping smotifs
        # ************************************************

        rmsd, transformed_coos = qcp.rmsdQCP(psmotif[0], csmotif_data[i], direction)

        if rmsd <= exp_data['rmsd_cutoff'][1]:
            # Loop constraint restricts the overlapping smotifs is not drifted far away.
            loop_constraint = llc.loopConstraint(transformed_coos, sse_ordered, direction, smotif_def)
            if loop_constraint:
                # Check whether the SSEs with in the assembled smotifs are clashing to one another
                no_clashes = qcp.clahses(transformed_coos, exp_data['clash_distance'])
            else:
                no_clashes = False


        if rmsd <= exp_data['rmsd_cutoff'][1] and no_clashes:
            # Prepare temp log array to save data at the end
            tlog, noe_fmeasure, total_percent, pcs_tensor_fits, rdc_tensor_fits = [], [], [], [], []

            tlog.append(['smotif', csmotif_data[i]])
            tlog.append(['smotif_def', sse_ordered])

            tlog.append(['qcp_rmsd', transformed_coos, sse_ordered, rmsd])
            cathcodes = sm.orderCATH(psmotif, csmotif_data[i][0], direction)
            tlog.append(['cathcodes', cathcodes])

            # ************************************************
            # Sequence filter
            # Aligns the smotif seq to target seq and calculates
            # sequence identity and the alignment score
            # ************************************************

            csse_seq, seq_identity, blosum62_score = Sfilter.S2SequenceSimilarity(current_ss, csmotif_data[i],
                                                                                  direction, exp_data)

            # concat current to previous seq
            concat_seq = sm.orderSeq(psmotif, csse_seq, direction)
            g_seq_identity = Sfilter.getGlobalSequenceIdentity(concat_seq, exp_data, sse_ordered)
            tlog.append(['seq_filter', concat_seq, csse_seq, seq_identity, blosum62_score])


            # ************************************************
            # Ambiguous NOE score filter
            # uses experimental ambiguous noe data to filter Smotifs
            # scoring based on f-measure?
            # ************************************************

            if 'noe_data' in exp_data_types:
                # noe_fmeasure = Nfilter.s3NOEfit(transformed_coos, sse_ordered, current_ss, exp_data)
                noe_fmeasure, no_of_noes = Noe.sXNOEfit(transformed_coos, sse_ordered, current_ss, exp_data)
                total_percent, total_noes = Noe.sXGlobalNOEfit(transformed_coos, sse_ordered, current_ss, exp_data)
                tlog.append(['NOE_filter', noe_fmeasure, no_of_noes])
                tlog.append(['GlobalNoe_filter', total_percent, total_noes])
            else:
                noe_fmeasure = False

            # ************************************************
            # Residual dipolar coupling filter
            # uses experimental RDC data to filter Smotifs
            # scoring based on normalised chisqr
            # ************************************************

            if 'rdc_data' in exp_data_types:
                rdc_tensor_fits = Rfilter.RDCAxRhFit2(transformed_coos, sse_ordered, exp_data, stage=2)
                if rdc_tensor_fits:
                    tlog.append(['RDC_filter', rdc_tensor_fits])
                else:
                    continue

            """
            if 'rdc_data' in exp_data_types:
                # if noe_fmeasure and noe_fmeasure >= exp_data['noe_fmeasure'][1]:
                if total_percent and total_percent >= exp_data['noe_fmeasure'][1]:
                    rdc_tensor_fits = Rfilter.RDCAxRhFit2(transformed_coos, sse_ordered, exp_data, stage=2)
                    tlog.append(['RDC_filter', rdc_tensor_fits])
                elif g_seq_identity >=100:
                    rdc_tensor_fits = Rfilter.RDCAxRhFit2(transformed_coos, sse_ordered, exp_data, stage=2)
                    tlog.append(['RDC_filter', rdc_tensor_fits])
                else:
                    continue
            """

            # ************************************************
            # Pseudocontact Shift filter
            # uses experimental PCS data to filter Smotifs
            # scoring based on normalised chisqr
            # ************************************************

            if 'pcs_data' in exp_data_types:
                pcs_tensor_fits = Pfilter.PCSAxRhFit2(transformed_coos, sse_ordered, exp_data, stage=2)
                tlog.append(['PCS_filter', pcs_tensor_fits])

            # ************************************************
            # Calc RMSD of the reference structure.
            # Used to identify the lowest possible RMSD
            # structure for the target, from the Smotif library.
            # ************************************************

            if 'reference_ca' in exp_data_types:
                ref_rmsd = ref.calcRefRMSD2(exp_data['reference_ca'], sse_ordered, transformed_coos, rmsd_cutoff=50.0)
            tlog.append(['Ref_RMSD', ref_rmsd, g_seq_identity])

            # if pcs_tensor_fits or rdc_tensor_fits:
            if pcs_tensor_fits:
                #dump data to the disk
                dump_log.append(tlog)

    # prevent dumping empty arrays with no data
    if len(dump_log) > 0:
        if 'rank_top_hits' in exp_data_types:
            dump_log = rank.rank_dump_log(dump_log, exp_data, stage=2)

        print "num of hits", len(dump_log),
        io.dumpPickle("tx_" + str(index_array[0]) + "_" + str(index_array[1]) + ".pickle", dump_log)

    return True
