#!/usr/bin/env python

# modification of Ali Fathi's code to parse Dimont's motif XML files (to work with mononucleotide PWMs instead of dinucleotide PWMs)

from bs4 import BeautifulSoup
import xml.etree.ElementTree as ET
import pandas as pd
import sys


MONONUCLEOTIDES = ['A', 'C', 'G', 'T']  # XY meand Y|X


def parse_parameters(param_tree):
    """
    children of param_tree:
    className, length, pos, pos, pos, pos (#length times)
    """
    param_tree_values = []
    for _pos in param_tree.contents[0:]:
        """
        className {}
        parameter {}
        [empty]
        """
        child_param = _pos[1]
        """
        [empty]
        value {}
        symbol {}
        index {}
        pseudoCount {}
        position {}
        context {}
        count {}
        free {}
        z {}
        t {}
        """
        # value = child_param.contents[1]  # className & "REAL SCORE!"
        # score = float(value.contents[1])
        score = float(child_param)
        param_tree_values.append(score)
    return param_tree_values


def read_xml(xml_file, output_folder, motif_name):
    with open(xml_file, 'r') as f:
        xml_motif = BeautifulSoup(f.read(), "lxml")
    # tree = ET.parse(xml_file)
    # root = tree.getroot()
    """
    alphabetcontainer --> className,
    length --> className,
    classWeights --> className {} - length {} - pos {'val': '0'} - pos {'val': '1'},
    params --> className {} - sequenceScoringParameterSet {},
    hasBeenOptimized --> className {},
    lastScore --> className {},
    score --> className {} - length {} - pos {'val': '0'} - pos {'val': '1'},
    beta --> className {} - length {} - pos {'val': '0'} - pos {'val': '1'} - pos {'val': '2'},
    prior --> class {}
    """
    # score = root[6]
    score = xml_motif.find_all("score")[0]
    """
    className {}
    length {} 2
    pos {'val': '0'}
    pos {'val': '1'}
    """
    # pos_1 = score[2]
    pos_1 = score.contents[2]
    """
    className {}
    ThresholdedStrandChIPper {}
    [empty]
    """
    # thresh = pos_1[1]
    thresh = pos_1.contents[1]
    """
    [empty]
    length {}
    starts {}
    freeParams {}
    function {}
    optimizeHidden {}
    plugIn {}
    hiddenParameter {}
    threshold {}
    """
    # func = thresh[3]
    func = thresh.contents[4]
    """
    className {}
    length {}
    pos {'val': '0'}
    """
    # pos = func[2]
    pos = func.contents[2]
    """
    className {}
    MarkovModelDiffSM {}
    [empty]
    """
    # markov = pos[1]
    markov = pos.contents[1]
    """
    [empty]
    bayesianNetworkSF {}
    [empty]
    lengthPenalty {}
    """
    # bayes = markov[0]
    bayes = markov.contents[1]
    """
    [empty]
    alphabets {}
    length {}
    trees {}
    isTrained {}
    ess {}
    numFreePars {}
    nums {}
    structureMeasure {}
    plugInParameters {}
    order {}
    roots {}
    freeParams {}
    """
    # trees = bayes[2]
    trees = bayes.contents[3]
    motif_length = int(trees.contents[1].text) - 1  # 15
    mononucleotide_df = pd.DataFrame(index=MONONUCLEOTIDES,
                                     columns=[str(i+1) for i in range(motif_length)])
    for position in range(motif_length+1):
        # className & parameterTree & [empty]
        param_pos = trees.contents[2+position]
        param_pos_tree = param_pos.contents[1]
        """
        [empty]
        pos {}
        contextPoss {}
        root {}
        firstParent {}
        firstChildren {}
        """
        param_pos_root = param_pos_tree.contents[3]  # className & treeElement & [empty]
        param_pos_tree = param_pos_root.contents[1]
        """
        [empty]
        contNum {}
        contextPos {}
        children {}
        pars {}
        """
        # if position == 0:  # pars
        # base_values = parse_parameters(
        #     param_tree=param_pos_tree.contents[4])
        # base_values_df = pd.DataFrame(base_values,
        #                               index=['A', 'C', 'G', 'T'],
        #                               columns=['0'])
        # else:  # children
        # className & treeElement
        if position > 0:
            param_pos_children = param_pos_tree.contents[4]
            mononuc_values = []
        #    Iterate over A, C, G, T
            for _child_pos in param_pos_children.contents[2:]:
                # <classname> & <treeelement> & [empty]
                _pos_tree = _child_pos.contents[1]
                """
                [empty]
                <contnum>
                <contextpos>
                <children>
                <pars>
                """
                _child_pos_pars = _pos_tree.contents[1]
                """
                <classname>
                <length>
                <pos val="0">
                <pos val="1">
                <pos val="2">
                <pos val="3">
                """
                # position_values = parse_parameters(param_tree=_child_pos_pars)
                position_values = float(_child_pos_pars.contents[1])
                mononuc_values.append(position_values)
            mononucleotide_df[str(position)] = mononuc_values
    if output_folder:
        if motif_name:
            # base_filename = motif_name + '_base_scores.csv'
            # mononucleotide_filename = motif_name + '_mononucleotide_scores.csv'
            mononucleotide_filename = motif_name + '.pwm'
        else:
            base_filename = 'base_scores.csv'
            mononucleotide_filename = 'mononucleotide_scores.csv'
        # base_values_df.to_csv(output_folder + '/' +
        #                       base_filename, sep='\t', index=True)
        mononucleotide_df.to_csv(output_folder + '/' +
                                 mononucleotide_filename, sep='\t', index=True)
    else:
        # print("Base values:")
        # print(base_values_df)
        print("mononucleotide PWM:")
        print(mononucleotide_df)


def main(args):
    # read_xml(xml_file="/Users/alifathi/Desktop/Codebook/Aim 1/CTCF_model.xml")
    read_xml(xml_file=args['input'],
             output_folder=args['output'], motif_name=args['name'])


def parse_args(args):
    if len(args) == 3 and args[1] == '-i':
        return {'input': args[2], 'output': None, 'name': None}
    elif len(args) == 5 and args[1] == '-i' and args[3] == '-o':
        return {'input': args[2], 'output': args[4], 'name': None}
    elif len(args) == 7 and args[1] == '-i' and args[3] == '-o' and args[5] == '-n':
        return {'input': args[2], 'output': args[4], 'name': args[6]}
    else:
        print('Usage: python DimontMotifParser.py -i <XML file> -o <output folder> -n <motif name>')
        sys.exit(1)


if __name__ == '__main__':
    args = sys.argv
    main(args=parse_args(args))