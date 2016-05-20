# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
import os
sys.path.append(os.path.join(os.environ['HOME'], 'dev/py'))
print sys.path

# <codecell>

from daltools import sirifc, one

# <codecell>

cmo = sirifc.sirifc(name="h2o/hf_h2o.SIRIFC").cmo

# <codecell>

S = one.read(label="AOONEINT", filename="h2o/hf_h2o.AOONEINT")
print S

sys.exit()

# <codecell>

print cmo.nrow, S.nrow

# <codecell>

print cmo.unblock()

# <codecell>

print S.unblock()

# <codecell>


