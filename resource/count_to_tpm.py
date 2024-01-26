import pandas as pd
import argparse as arg
if __name__ == "__main__":
    parser = arg.ArgumentParser()
    parser.add_argument('-f', '--file', type=str, action="store",default=False)
    args = parser.parse_args()
    with open(args.file) as f:
        first_line = f.readline()
    count_df = pd.read_csv(args.file,sep = "\t",skiprows=1)
    count_df['TPM'] = ((count_df.iloc[:,6]*1000)/count_df.iloc[:,5])/(((count_df.iloc[:,6]*1000)/count_df.iloc[:,5]).sum())*1000000
    count_df.to_csv(args.file,sep = "\t",index=0)
    with open(args.file, 'r+') as f:
        content = f.read()        
        f.seek(0, 0)
        f.write(first_line+content)

