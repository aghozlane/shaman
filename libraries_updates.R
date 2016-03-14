#### to avoid conflict between d3heatmap and rnvd3 #####
# the conflict come from the fact that both libraries use d3 but not the same version
# when the most recent one is loaded (the one of d3heatmap), it breaks the one of rNVD3
# the idea here is to modify the d3 of d3heatmap so that the library create an object called d33
# then modify the d3heatmap files so that it use this d33 object instead of d3

# first get d3heatmap javascript library path
d3heatmap_path = paste(.libPaths(), 'd3heatmap/htmlwidgets/lib/', sep='/')

# then get all js files in this directory
js_files = dir(path=d3heatmap_path, pattern = "*.js$", recursive = T, full.names = T)

# then create a function tu update the files
d3_to_d33 = function(target){
  # the function will use the bash tool sed on the file target
  # it makes 3 replacements
  system(command=paste('sed -i .bck -e "s/d3\\ /d33 /g" -e "s/d3\\./d33./g" -e "s/d3;/d33;/g" ', target))
}
# apply the replacement function to all js files
sapply(js_files, d3_to_d33)

# change from d3.min.js to d3.js in the yaml file
system(command=paste('sed -i .bck -e "s/d3\\.min/d3/g" ', .libPaths(), '/d3heatmap/htmlwidgets/d3heatmap.yaml', sep=''))

# it's done. d3heatmap should use d33 object now

#### replacement of the css cals of the rNVD3 ####

# get the rNVD3 path
rNVD3_path = paste(.libPaths(), 'rNVD3/nvd3/css/', sep='/')

# get the path of the current source. may need to be modified (parent.frame(1) or parent.frame(2) or parent.frame(3) or else) 
current_path = dirname(parent.frame(2)$ofile) 
# get the css files path
css_files = dir(path=current_path, pattern = "*.css$", full.names = T)
# create the function to copy the css file
cp_to_rNVD3 = function(p, dest){
  system(paste('cp', p, dest))
}
#apply the copy to all the css files
sapply(css_files, cp_to_rNVD3, rNVD3_path)

# first get d3heatmap javascript library path
scatterD3_path = paste(.libPaths(), 'scatterD3/htmlwidgets/', sep='/')

# then get all js files in this directory
js_files = dir(path=scatterD3_path, pattern = "*.js$", recursive = T, full.names = T)

# then create a function tu update the files
d3_to_d333 = function(target){
  # the function will use the bash tool sed on the file target
  # it makes 3 replacements
  system(command=paste('sed -i .bck -e "s/d3\\ /d333 /g" -e "s/d3\\./d333./g" -e "s/d3;/d333;/g" -e "s/d3=/d333=/g" ', target))
}
# apply the replacement function to all js files
sapply(js_files, d3_to_d333)

# change from d3.min.js to d3.js in the yaml file
system(command=paste('sed -i .bck -e "s/d3$/d333/g" ', .libPaths(), '/scatterD3/htmlwidgets/scatterD3.yaml', sep=''))

