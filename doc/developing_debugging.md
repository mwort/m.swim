# Using the central m.swim.git repository

## Requirements and copying the m.swim code

Frist of all you need to have git installed on your machine or you can use the cluster git.
With a direct cluster connection, you can copy the m.swim code to your machine like this:

    git clone git@gitlab.pik-potsdam.de:wortmann/m.swim.git


## Working with git and changing the code

Now you have the entire code (incl. documentation) in the m.swim directory wherever you executed the clone command. Make your changes to the code and commit them inside local m.swim directory like this:
Check what files are changed and will be commited

    git status

Check what lines have changed

    git diff

Show what has been prevously changed (the commits) and branches (opens window with gui, i.e. you need X Windows set up!)

    gitk

## Committing changes locally

Add individual files to the commit

    git add <file>

... and then commit

    git commit -m 'Short description of what has changed'

... or you just commit all files that have changed without adding them:

    git commit -a -m 'Short description of what has changed'

## Committing changes to the cluster m.swim repository

Now you can make these commits to the central cluster directory by 'pushing' them to it. You'll first have to set up the remote repository with any name (pik.m.swim) like this:

    git remote add pik.m.swim <username>@cluster.pik-potsdam.de:/home/wortmann/src/grass7/m.swim.git

Or you can leave out the '<username>@cluster.pik-potsdam.de:' if you are working on the cluster.
Now you can push the master branch (default, could be any other) to the pik.m.swim repository:

    git push pik.m.swim master

## Updating your local m.swim code from the central m.swim repository

To update your local m.swim repository just 'fetch' or 'pull' (fetch+merge) it using the previously set up server name pik.m.swim:
To first fetch the changes and examine them

    git fetch pik.m.sim master

Or to fetch and merge the changes to your local repository straight away

    git pull pik.m.swim master
