

--- GIT instructions:

When you make changes to this folder, sync those changes with the repository in GitHub. This is how to do it:

# Open Git Bash
# Change directory to this folder by writing:
	cd /e/Dropbox/NMDA
# You can see the contents of the folder by writing:
	ls
# Stage all the changes you want to save:
	git add .
# Commit changes to your local repository (this is similar to the traditional "Save" function)
	git commit -m "Message explaining the changes you made since the last saved version"
# Check the status of what you've done or the list of commits by writing, respectively:
	git status
	git log
# Now, you want to push those changes to the repo on GitHub:
	git push origin master