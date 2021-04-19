"""Build LAMMPS input file
..todo:: add error checking    
"""
from mdlib.pre import build
from mdlib.post import fixphonon
from argparse import ArgumentParser
    
if __name__=="__main__":

    parser = ArgumentParser(description="Helper module for the dispersion curve md simulation")
    parser.add_argument("mode", type=str, choices=["pre", "post"])
    parser.add_argument("collect", nargs="+")
    args = parser.parse_args()

    if args.mode == "pre":
        print("... pre-processing")
        fin = str(args.collect[0])
        nx = int(args.collect[1])
        ny = int(args.collect[2])
        nz = int(args.collect[3])
        fout = str(args.collect[4])
        print("    ", fin,nx,ny,nz,fout)
        build.write_from_file(fin,nx,ny,nz,fout)

    if args.mode == "post":
        print("... post-processing")
        wdir = str(args.collect[0])
        pattern = str(args.collect[1])
        print("    ", wdir, pattern)
        files, pvals = fixphonon.find_files(wdir,pattern)
        print("    ",len(files)," files")
        qr, ddisp, dintp, labels = fixphonon.process_log_files(files)
        fixphonon.make_interactive_plot(qr,pvals,ddisp,dintp,labels)
