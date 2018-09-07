#########################################################################
# File Name: install.sh
# Created Time: Thu 06 Sep 2018 07:08:08 PM CST
#########################################################################
#!/bin/bash
ln -s src/Basic.py basic.py



# install somatic
ln -s ../src/somatic/pipe.py bin/somatic.py
chmod 755 src/somatic/maftag.py 
chmod 755 src/somatic/vcfaddtag.py bin/somatic.py
chmod 755 src/somatic/vardict_msi2.py
