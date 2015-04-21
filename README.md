#  mirnapy  
### Author: Gene Blanchard  
### Email: me@geneblanchard.com  

## Requirements  
*  fastx_toolkit  
*  bwa   

If using *ubuntu* you can install with `sudo apt-get install bwa fastx-toolkit`, you may also need to install `libgtextutils`. 

If using  *mac* I recomend recomend to install [homebrew](http://brew.sh/). Next install the *homebrew-science* repo with `brew tap homebrew/science`. Install the *fastx-toolkit* and *bwa* with `brew install bwa fastx_toolkit`

## Example 
`mirna_pipe.py -i A1.fastq,B1.fastq,A2.fastq,B2.fastq -o A1_B1_A2_B2_output`  
This would create the `A1_B1_A2_B2_output` folder in your current working directory.  

# Database
The default database is located at `DB_mature/mouse/mature_dna_mouse.fa` change this with the `-d` option

## Output 
The output contains the folowing:
*  **commands.txt**  
  *  A file that lists all commands run  
*  **clipped_files/**  
  *  The output of the fastx clipping   
*  **alignment_files/**  
  *  The output of the BWA alignment  
*  **SEQEM/**  
  *  The results of the SEQEM command  
*  **counts/**  
  *  Raw counts from the SEQEM command  
*  **RSEM/**  
  *  Groomed tab-seperated files that are ready to input into RSEM  


# TODO
Need to add sequence length filtering options
