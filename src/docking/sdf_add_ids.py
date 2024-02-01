import argparse

def get_ids(sdf):

    id_header = False
    ids = []

    with open(sdf) as f:
        for line in f:
            if id_header:
                ids.append(line.strip())
                id_header = False
            if line.startswith("> <chembl_id>"): id_header = True

    return ids

def insert_ids(sdf_in, ids, sdf_out=None):
    
    if not sdf_out: sdf_out = sdf_in.split(".")[0]+"_with_ids.sdf"
    block_start=True
    file_out = open(sdf_out,"w")
    with open(sdf_in,"r") as f:
        for line in f:
            if block_start: 
                file_out.write(ids.pop(0)+"\n")
                block_start = False
            else:
                file_out.write(line)
                if line.startswith("$$$$"): block_start = True
    file_out.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sdf", type=str, help="sdf file from chembl download")
    parser.add_argument("-o", "--out-file", dest="out_file", help="output sdf (with ids) name")
    args = parser.parse_args()

    ids = get_ids(args.sdf)
    insert_ids(args.sdf, ids, args.out_file)

if __name__=="__main__":
    main()
