# The functions here were copied from the Xyce repository:
# https://github.com/Xyce/Xyce/tree/master/utils

import os, sys, re
import numpy

#-------------------------------------------------------------------------------
def findBlock(text, pos=0, delim="'"+'"'+"\(\[\{"):
    """Given the input text (potentially multiline) and an optional pos marking
    the starting position, find an opening delimeter -- either (, [, {, single
    quote, or double quote -- and return a tuple of integers indicating the
    character indexes of the text block -- closed with ), ], }, single quote, or
    double quote, respectively -- while correctly handling nested blocks."""

    # Define delimeter strings based on optional delim input
    quote1Delimeter = ""
    quote2Delimeter = ""
    if "'" in delim:
      quote1Delimeter += "'"
    if '"' in delim:
      quote2Delimeter += '"'
    openDelimeters = ""
    closeDelimeters = ""
    if "\(" in delim:
      openDelimeters  += "\("
      closeDelimeters += "\)"
    elif "\[" in delim:
      openDelimeters  += "\["
      closeDelimeters += "\]"
    elif "\{" in delim:
      openDelimeters  += "\{"
      closeDelimeters += "\}"

    import re
    # Define delimeter regular expressions
    if quote1Delimeter:
      quote1RE = re.compile("([" + quote1Delimeter + "])", re.M)
    if quote2Delimeter:
      quote2RE = re.compile("([" + quote2Delimeter + "])", re.M)
    openRE   = re.compile("([" + openDelimeters  +
                                 quote1Delimeter +
                                 quote2Delimeter + "])", re.M)
    anyRE    = re.compile("([" + openDelimeters  +
                                 quote1Delimeter +
                                 quote2Delimeter +
                                 closeDelimeters + "])", re.M)

    # Find the first opening delimeter
    matchObject = openRE.search(text, pos)
    if not matchObject: return (None, None)

    # Initialize the loop
    stack = [ matchObject.group() ]
    start = matchObject.start()
    pos   = start + 1

    # Find the end of the block
    while matchObject:

        # Determine the active delimeter regular expression
        if   stack[-1] == quote1Delimeter:
            activeRE = quote1RE
        elif stack[-1] == quote2Delimeter:
            activeRE = quote2RE
        else:
            activeRE = anyRE

        # Search for the next delimeter
        matchObject = activeRE.search(text, pos)
        if matchObject:
            delimeter = matchObject.group()
            pos       = matchObject.end()

            # Check for matched delimeters
            if (((stack[-1] == quote1Delimeter) and
                 (delimeter == quote1Delimeter)) or
                ((stack[-1] == quote2Delimeter) and
                 (delimeter == quote2Delimeter)) or
                ((stack[-1] == "("            ) and
                 (delimeter == ")"            )) or
                ((stack[-1] == "["            ) and
                 (delimeter == "]"            )) or
                ((stack[-1] == "{"            ) and
                 (delimeter == "}"            ))   ):
                stack.pop()                  # Remove the last element from the list
                if len(stack) == 0:
                    return (start, pos)

            # Process unmatched delimeter
            else:
                if (delimeter in openDelimeters  or
                    delimeter == quote1Delimeter or
                    delimeter == quote2Delimeter   ):
                    stack.append(delimeter)  # Add the delimeter to the stack
                else:
                    raise RuntimeError("findBlock: mismatched delimeters: " + \
                          stack[-1] + " " + delimeter)

    # We made it through all of text without finding the end of the block
    raise RuntimeError("findBlock: open block: " + join(stack))

#-------------------------------------------------------------------------------
def getXyceData(file,verbose=False):
  """
  getdata(file) reads the data from a Xyce output prn file into a list of
  column names and an array containing the data.
  """
  if os.path.exists(file):
    input = open(file,'r').readlines()
    numlines = len(input)-2
    # remove spaces and braces from expressions
    line0 = input[0]
    match = findBlock(line0,delim="\{")
    while match[0]:
      beg = match[0]
      end = match[1]
      group = line0[beg+1:end-1] # remove braces
      line0 = line0[:beg] + re.sub(r"[ ]",r"",group) + line0[end+1:]
      match = findBlock(line0,end,delim="\{")
    tags = line0.split()
    numentries = len(tags)
    if verbose:
      print ("Read file: " + file)
      print ("Found columns:  ", end=' ')
      print (tags)
      print ("Found " + str(numlines) + " lines of data" )
    data = numpy.zeros((numlines,numentries),'double')
    if verbose:
      print ("Reading lines: ",)
    for i in range(1,len(input)-1):
      if verbose:
        print (".", end=' ')
      s = input[i].split()
      if s[0] != 'End':  # if this is the final line text string, skip it.
        if len(s) != numentries:
          print ("Error: " + file + ":" + str(i+1))
          print ("Number of columns read is not equal to number of columns in header.")
          print ("numentries = " + str(numentries))
          print ("len(s) = " + str(len(s)))
          sys.exit(1)
        data[i-1,:] = [float(j) for j in s]
    if verbose:
      print ("\n")
  else:
    print ("Error, file, " + file + " does not exist")
    print (getXyceData.__doc__)
    sys.exit(1)
  return (tags,data)