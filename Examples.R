
library(mmp)

library(panelr)
library(rms)
library(lmtest)
library(multiwayvcov)
library(lfe)
data("WageData")


mod1 <- lm(lwage~exp+fem, WageData)
mod2 <- lm(lwage~exp+fem+occ+ind, WageData)
mod3 <- lm(lwage~exp+fem+occ+ind+ms, WageData)
WageData$id <- factor(WageData$id)
mod4 <- lm(lwage~exp+fem+occ+ind+ms+id, WageData)
mod5 <- lm(lwage~exp*fem+occ+ind+ms+id, WageData)
WageData$IncomeGroups <- factor((WageData$lwage < median(WageData$lwage))*1, c(0, 1), c("Low Income", "High Income"))
mod67 <- lapply(split(WageData, WageData$IncomeGroups), function(x)lm(lwage~exp+fem, x))
mod8 <- lrm(fem~lwage+exp+t, WageData, x = T, y= T)
mod8 <- robcov(mod8, WageData$id)
mod9 <- lm(lwage~exp+fem+occ+ind+ms+blk, WageData)
mod9c <- coeftest(mod9, cluster.vcov(mod9, cluster = WageData$id))
mod10 <- felm(lwage~exp+fem|t|0|id, WageData)



# Example 1
View(mmp(list(mod1, mod2, mod3)))

# Example 2
View(mmp(list(Baseline = mod1, "Some Controls" = mod2, "More Controls" = mod3, "Awful" = mod4)))
# Example 3
View(mmp(list(Baseline = mod1, "Some Controls" = mod2, "More Controls" = mod3, "Fixed Effects" = mod4, "Interaction" = mod5),
         CoefNames = list("Female first" = "fem", "<Generic>", "Experience at the end" = "exp"),
         Factors = list("Individual FE" = "id"),
         BottomMatter = c("N", "r.squared", "f-statistics", "DF")))
# Example 4
View(mmp(mod67))
# Example 5
View(mmp(list(M8 = mod8, M9 = mod9, M10 = mod10),
         CoefNames = list("<Generic>", "Uninteresting variables" = "(ind)|(ms)"),
         CoefTests = list(NULL, mod9c), ShowDependentVariable = T, ShowModelType = T,
         BottomMatter = c("StandardErrorInfo")))

View(mmp(list(M8 = mod8, M9 = mod9, M10 = mod10),
         CoefNames = list("<Generic>", "Uninteresting variables" = "(ind)|(ms)"),
         CoefTests = list(NULL, mod9c), ShowDependentVariable = T, ShowModelType = T,
         Factors   = list("Uninteresting variables"),
         BottomMatter = c("StandardErrorInfo")))

# Example 6
# write.table(mmp(list(M8 = mod8, M9 = mod9, M10 = mod10),
#               CoefNames = list("<Generic>", "Uninteresting variables" = "(ind)|(ms)"),
#               CoefTests = list(NULL, mod9c), ShowDependentVariable = T, ShowModelType = T,
#               Factors   = list("Uninteresting variables"),
#               BottomMatter = c("StandardErrorInfo")), file =  "test.csv", sep = ";")

