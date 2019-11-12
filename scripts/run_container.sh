#!/usr/bin/env bash
set -e

# adapted from: http://wiki.ros.org/docker/Tutorials/GUI

# Find the invoking user's group ID and set the group ID of the user inside the
# container to the same value. There's a very good chance that this will result
# in an error complaining that the container cannot find the name of this group,
# but that's OK.
IMG_USR=root
IMG_GID=$(id -g)

if [ $# != 0 ]
then
    task="$@"
else
    task=/bin/bash
fi

docker run \
    -it \
    --net=host \
    --user=$IMG_USR:$IMG_GID \
    --env="DISPLAY" \
    --workdir="/home/$IMG_USR" \
    --volume="/home/$USER/Develop:/Develop" \
    --volume="/etc/shadow:/etc/shadow:ro" \
    --volume="/mnt/prostor1:/mnt/prostor1:rw" \
    --volume="/mnt/cds-scratch/:/mnt/cds-scratch/:rw" \
    cmriprocan/r-base \
    /bin/bash
