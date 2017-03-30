import sys
import pstats

p = pstats.Stats(sys.argv[1])

p.sort_stats('filename')
p.print_stats()


