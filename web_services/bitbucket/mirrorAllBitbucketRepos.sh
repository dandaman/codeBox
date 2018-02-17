#!/bin/bash
#Script to get all repositories under a user from bitbucket
#Usage: getAllRepos.sh login user_or_team_name
#source://movingtothedarkside.wordpress.com/2015/01/10/clone-all-repositories-from-a-user-bitbucket/
#modified by DL
#17.02.2018

curl -u ${1} https://api.bitbucket.org/1.0/users/${2} > repoinfo

for repo_name in `cat repoinfo | sed -r 's/("name": )/\n\1/g' | sed -r 's/"name": "(.*)"/\1/' | sed -e 's/{//' | cut -f1 -d\" | tr '\n' ' '`
do
    echo "Cloning " $repo_name
    git clone --mirror ssh://git@bitbucket.org/${2}/$repo_name
    echo "---"
done
