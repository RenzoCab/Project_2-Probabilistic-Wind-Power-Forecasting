import argparse
parser = argparse.ArgumentParser(description='Parsing input configuration file...')
parser.add_argument('-filename', help=' Config file name or path')
parser.add_argument('--version', action='version', version='Likelihood Evaluato 1.0')

args = parser.parse_args()
print(args.filename)
