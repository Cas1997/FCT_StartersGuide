#!/bin/bash

export ALIBUILD_WORK_DIR=/misc/alidata150/alice_u/cas/O2Simulation/alice/sw
eval `/usr/local/bin/alienv shell-helper`
export PATH=${PATH}

# o2_env.sh
echo "#!/bin/bash" > o2_env.sh
echo "" >> o2_env.sh
alienv printenv O2/latest-retractable_IRIS-o2 >> o2_env.sh
sed -i "s|test 0;||g" o2_env.sh
chmod a+x o2_env.sh
