import os,glob,string,numpy

def read_file2(filename,delimiter='\t'):
    """ Alternativ: numpy.loadtxt(filename,dtype=numpy.int64) bzw.
        dt={'names':('time','infected'),'formats':(np.int64,np.float32)}
        dann: numpy.loadtxt(filename,dtype=dt)
    """
    
    # Liest die Datei ein
    file_tab=file(filename,"r")
    line=file_tab.readline()
    data=[]
    while line!="":
        row=[]
        for val in string.split(line[:-1],delimiter):
            if string.strip(val)!="":
                row.append(string.strip(val))
        data.append(row)    
        line=file_tab.readline()
    file_tab.close()            
    return data

def merge_txt_files(outfile='Merged.dat',subdir=""):
    """ Merge multiple txt-files in working directory.
        Subdir: "/my_subdirectory"
        
        Returns joint txt-File
    """
    if subdir=="" and outfile[-4:] == ".txt":
        print 'Will convert ',outfile,' to ',outfile[0:-4],'.dat'
        outfile=outfile[0:-4]
        outfile+='.dat'
    
    merged=[]
    for file in glob.glob(os.path.join(subdir+"*.txt")):
        data=read_file2(file) #numpy.loadtxt(file,dtype=datatype)
        merged.extend(data)
    
    numpy.savetxt(outfile,merged,delimiter='\t',fmt='%s')
    
    return

if __name__=="__main__":
    merge_txt_files()
    
