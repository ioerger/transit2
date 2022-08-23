methods = {}
from .combined_wig import Export; methods["combined_wig"] = Export()
from .igv          import Export; methods["igv"         ] = Export()
from .mean_counts  import Export; methods["mean_counts" ] = Export()
# from .norm         import Export; methods["norm"        ] = Export()