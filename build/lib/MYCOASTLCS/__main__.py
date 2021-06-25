# -*- coding: utf-8 -*-

from .Common import Common
import argparse

def printMyCoast():
    print("""    
     __          __    _                                _                              
     \ \        / /   | |                              | |                             
      \ \  /\  / /___ | |  ___  ___   _ __ ___    ___  | |_  ___                       
       \ \/  \/ // _ \| | / __|/ _ \ | '_ ` _ \  / _ \ | __|/ _ \                      
        \  /\  /|  __/| || (__| (_) || | | | | ||  __/ | |_| (_) |                     
         \/  \/  \___||_| \___|\___/ |_| |_| |_| \___|  \__|\___/                      
                                                                                       
                                                                                       
      __  __          _____                    _    _       _____   _____    ___    __ 
     |  \/  |        / ____|                  | |  | |     / ____| / ____|  / _ \  /_ |
     | \  / | _   _ | |      ___    __ _  ___ | |_ | |    | |     | (___   | | | |  | |
     | |\/| || | | || |     / _ \  / _` |/ __|| __|| |    | |      \___ \  | | | |  | |
     | |  | || |_| || |____| (_) || (_| |\__ \| |_ | |____| |____  ____) | | |_| |_ | |
     |_|  |_| \__, | \_____|\___/  \__,_||___/ \__||______|\_____||_____/   \___/(_)|_|
               __/ |                                                                   
              |___/                                                                    
    """)
    
    print('Authors')
    print('*Angel Daniel Garaboa Paz. Group of Non Linear Physics')
    print('*Vicente Perez Muñuzuri. Group of Non Linear Physics')
    print('*James Clark. Plymouth Marine Laboratory')
    print('*Ricardo Torres. Plymouth Marine Laboratory')
    print('\n')
                                                        

def main():
    
    printMyCoast()
    
    argParser = argparse.ArgumentParser(description='Indexes input files for MOHID Lagrangian to parse. Use -h for help.')
    argParser.add_argument("-j", "--json", dest="case_json",
                    help=".json file with the case definition for the MOHID Lagrangian run", metavar=".json")
    argParser.add_argument("-i", "--input", dest="input_file",
                    help="input netcdf file or ensemble of netcdf files", metavar="input_file")
    argParser.add_argument("-o", "--output", dest="output_file",
                    help="output netcdf file with desire fields", metavar="output_file")
    args = argParser.parse_args()
    
    run = Common()
    run.read_json(args.case_json)
    run.run_ftle_lcs(args.input_file,args.output_file)
    
    print('Finish!\n\n')

main()