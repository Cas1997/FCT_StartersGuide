#!/bin/bash

export ALIBUILD_WORK_DIR=<path to alibuild work dir>
eval `/usr/local/bin/alienv shell-helper`
export PATH=${PATH}

# o2_env.sh
echo "#!/bin/bash" > o2_env.sh
echo "" >> o2_env.sh
alienv printenv <your O2 environment> >> o2_env.sh
sed -i "s|test 0;||g" o2_env.sh
chmod a+x o2_env.sh
