
import sys

if (len(sys.argv) != 2):
  sys.stderr.write('\nRun the program: parintPAFtargets [PAF file].\n')
  exit();

paffilename = sys.argv[1]
targets_dict = {}

with open(paffilename, 'r') as paffile:
    for line in paffile:
      if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
        pass
      else:
        try:
            elements = line.split('\t')    # splitting with tab as delimitters
            elcount = len(elements)

            pafline = {}
            attributes = {}

            pafline['QNAME'] = elements[0]
            pafline['QLEN'] = int(elements[1])
            pafline['QSTART'] = int(elements[2])
            pafline['QEND'] = int(elements[3])
            pafline['STRAND'] = elements[4]
            pafline['TNAME'] = elements[5]
            pafline['TLEN'] = int(elements[6])
            pafline['TSTART'] = int(elements[7])
            pafline['TEND'] = int(elements[8])
            pafline['NRM'] = int(elements[9])
            pafline['ABL'] = int(elements[10])
            pafline['MQUAL'] = int(elements[11])

            if elcount > 12:
              for i in range(12, elcount):
                element = elements[i]
                att_list = element.split(':')
                attributes['TAG'] = att_list[0]
                attributes['TYPE'] = att_list[1]
                attributes['VALUE'] = att_list[2]
                pafline['ATTRIB'] = attributes

            targets_dict[pafline['TNAME']] = 1;
        except:
            import pdb
            pdb.set_trace()
            pass

for tname in targets_dict.iterkeys():
  sys.stdout.write(tname + '\n')
