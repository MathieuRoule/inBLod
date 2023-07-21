###########################
# Check for packages installation
###########################

using Pkg # For packages management

function pkgcheck(pkgnames...)

    #####
    # Listing all installed packages
    #####
    # See function Pkg.installed() (deprecated)
    installedpkg = Array{String,1}()
    deps = Pkg.dependencies()
    for (uuid, dep) in deps
        dep.is_direct_dep || continue
        append!(installedpkg,[dep.name])
    end

    #####
    # Listing the uninstalled needed packages
    #####
    unknownpkg = Array{String,1}()
    for pkg in pkgnames
        if !(pkg in installedpkg)
            append!(unknownpkg,[pkg])
        end
    end

    #####
    # Installing the unknown packages
    #####
    if !(isempty(unknownpkg))
        println("The following packages are needed:")
        for pkg in unknownpkg
            println(pkg)
        end
        print("Do you want to install these packages ([y]/n) ?")
        answer = chomp(readline())
        if (answer == "") ||(answer == "y") || (answer == "yes")
            for pkg in unknownpkg
                Pkg.add(pkg)
            end
        end
    end
end


