import damaris

f = file("/tmp/script_out.txt","a")

for v in damaris.get_variables():
  print >> f, "has var",v.name,
  chunks = v.chunks
  if len(chunks) == 0:
    print >> f, "no chunks"
  else:
    print >> f, len(chunks),"chunks"

  for c in v.chunks:
    print >> f, "it ",c.iteration, "from rank",c.source, "    dims",c.lower_bounds,c.upper_bounds,"= \n" #,c.data


print >> f, "Cleaning everything in the variables\n"
damaris.clear()

f.close()
