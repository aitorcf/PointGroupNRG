using PkgTemplates 

t = Template(;
    user = "aitorcf",
    authors = "Aitor Calvo-Fernández",
    plugins = [
               License(name="MIT"),
               Git(),
               GitHubActions(),
              ],
   )

t("PointGroupNRG")

