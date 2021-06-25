# -*- coding: utf-8 -*-

from .Common import Common
import argparse


MY_COAST_LOGO = """
███╗░░░███╗██╗░░░██╗░█████╗░░█████╗░░█████╗░░██████╗████████╗██╗░░░░░░█████╗░░██████╗
████╗░████║╚██╗░██╔╝██╔══██╗██╔══██╗██╔══██╗██╔════╝╚══██╔══╝██║░░░░░██╔══██╗██╔════╝
██╔████╔██║░╚████╔╝░██║░░╚═╝██║░░██║███████║╚█████╗░░░░██║░░░██║░░░░░██║░░╚═╝╚█████╗░
██║╚██╔╝██║░░╚██╔╝░░██║░░██╗██║░░██║██╔══██║░╚═══██╗░░░██║░░░██║░░░░░██║░░██╗░╚═══██╗
██║░╚═╝░██║░░░██║░░░╚█████╔╝╚█████╔╝██║░░██║██████╔╝░░░██║░░░███████╗╚█████╔╝██████╔╝
╚═╝░░░░░╚═╝░░░╚═╝░░░░╚════╝░░╚════╝░╚═╝░░╚═╝╚═════╝░░░░╚═╝░░░╚══════╝░╚════╝░╚═════╝░

██╗░░░██╗░█████╗░░░░░░███╗░░
██║░░░██║██╔══██╗░░░░████║░░
╚██╗░██╔╝██║░░██║░░░██╔██║░░
░╚████╔╝░██║░░██║░░░╚═╝██║░░
░░╚██╔╝░░╚█████╔╝██╗███████╗
░░░╚═╝░░░░╚════╝░╚═╝╚══════╝
"""


def printMyCoast():

    print(MY_COAST_LOGO)

    print('-> Authors:')
    print('-> Angel Daniel Garaboa Paz. Group of Non Linear Physics (GFNL).')
    print('-> Vicente Perez Muñuzuri. Group of Non Linear Physics (GFNL).')
    print('-> James Clark. Plymouth Marine Laboratory (PML).')
    print('-> Ricardo Torres. Plymouth Marine Laboratory (PML).')
    print('\n')


def main():

    printMyCoast()

    argParser = argparse.ArgumentParser(
        description='FTLE and LCS calculator for ensembles of Lagrangian \
            simulations. Use -h for help.')
    argParser.add_argument("-j", "--json",
                           dest="case_json",
                           help=".json file with setup fpr FTLE/LCS",
                           metavar=".json")
    argParser.add_argument("-i", "--input",
                           dest="input_file",
                           help="input netcdf file/s files",
                           metavar="input_file")
    argParser.add_argument("-o", "--output",
                           dest="output_file",
                           help="output netcdf file with desire fields",
                           metavar="output_file")
    args = argParser.parse_args()

    run = Common()
    run.read_json(args.case_json)
    run.run_ftle_lcs(args.input_file, args.output_file)

    print('Finish!\n\n')


main()
