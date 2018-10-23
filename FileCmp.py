import sys 


numFiles = len(sys.argv) - 1
files = (sys.argv)
files.pop(0) #remove the 1st element as it is the script name
print ('Number of files:', numFiles)
print ('Files names:', files)

same = True

if numFiles == 2:
    with open(files[0]) as f1, open(files[1]) as f2:
        f1Lines = sum(1 for line in f1)
        f2Lines = sum(1 for ling in f2)
        if (f1Lines == f2Lines):
            for l1,l2 in zip(f1,f2): 
                if (l1 != l2):
                    same = False
                    print (l1, l2)
                    break
        else: 
            print ("Files are not the same length")
            same = False

elif numFiles == 3:
    with open(files[0]) as f1, open(files[1]) as f2, open(files[2]) as f3:
        f1Lines = sum(1 for line in f1)
        f2Lines = sum(1 for ling in f2)
        f3Lines = sum(1 for ling in f3)
        if (f1Lines == f2Lines == f3Lines):
            for l1,l2,l3 in zip(f1,f2,f3):
                if (l1 != l2 or l1 != l3 or l2 != l3):
                    same = False
                    print (l1, l2, l3)
                    break
        else:
            print ("Files are not the same length")
            same = False 

if (same):
    print("The files are the same")
else:
    print("The files are not the same")

    