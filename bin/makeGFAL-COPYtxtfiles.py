# MakeGFAL-COPYTxtFiles 

# Get grid proxy first
#voms-proxy-init --voms cms

## Use this python file to make any gfal-copy file you want
# Options
# File location
# Txt File name
import argparse , os

# set up an argument parser                                                                                                                                                                         
parser = argparse.ArgumentParser()

parser.add_argument('--loc', dest = 'LOC', default= 'srm://maite.iihe.ac.be:8443/pnfs/iihe/cms/store/user/schenara/MC_RunII_2016/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/crab_DYJetsToLL_M-10to50_v3/191130_081246/0000/' )
parser.add_argument('--name', dest='NAME', default='2016_DYM10to50_ext0')# used to name txt file and eos directory 
parser.add_argument('--dest', dest='EOSDEST', default='/eos/user/a/asparker/TopLFV/')# eos directory to move files to
parser.add_argument('--v', dest='VERBOSE', default=True)
parser.add_argument('--cp', dest='CPTOEOS', default=True) # Move these files to my EOS space

ARGS = parser.parse_args()


# Check that files exist and you can contact their server location using gfal-ls
print "Test gfal-ls"

testlstxt = "gfal-ls " + ARGS.LOC

os.system( testlstxt )



# Make a text file containing the gfal-ls information
#myfile = open( ARGS.NAME + ".txt", 'w')

# Loop over gfal-ls output and add ARGS.LOC to file name so they are in correct format for gfal-copy
#gfal-copy -f --from-file files.txt file://$PWD 
# where files.txt is a file where every line is a source in srm url syntax


# create list of files in SRM URL syntax
txtfileList = []

for fn in filter(None,os.popen( testlstxt  ).read().split('\n')):
    #print "FILE X : {}".format(fn)
    txtfileList.append( ARGS.LOC + fn )
    if ARGS.VERBOSE :  print "Adding to txt file :  {}".format( ARGS.LOC + fn  ) 

if ARGS.VERBOSE :
    print "txtfileList"
    print txtfileList
 

#'\n'

# Write the list to a text file
with open( ARGS.NAME + ".txt", "w") as output:
    for f in txtfileList  :
        output.write(str(f) + '\n')



print"gfal-ls output written to file {}".format( ARGS.NAME + '.txt'  )
#myfile.close()




if ARGS.CPTOEOS :
    # Move these files to my EOS space
    #gfal-copy -f --from-file files.txt file://$PWD                                                                                             # where files.txt is a file where every line is a source in srm url syntax                                       
    isreal = os.system( "ls "+  ARGS.EOSDEST+ARGS.NAME  )  
    print  ARGS.EOSDEST+ARGS.NAME 
    print isreal
    
    if isreal > 0 :  
        print "create eos dir if it doesn't exist: {}".format(ARGS.EOSDEST+ARGS.NAME)
        os.system( "mkdir "+ ARGS.EOSDEST+ARGS.NAME )

        # Copy files to their new home
        gfcom = "gfal-copy -f --from-file  " + ARGS.NAME+ ".txt  " + ARGS.EOSDEST+ARGS.NAME + "/" 
        print gfcom
        os.system( gfcom )
        print "Done moving all files from {} to {}".format ( ARGS.LOC , ARGS.EOSDEST+ARGS.NAME  )
    elif isreal == 0 :
        print "eos dir exists, skipping..."
    
