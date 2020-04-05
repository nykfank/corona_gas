# COVID-19 pandemic visualisation by Niklaus Fankhauser
# git clone https://github.com/CSSEGISandData/COVID-19
# Requires: sf, ggplot2, rnaturalearth, rnaturalearthhires, rgeos, git2r

library("ggplot2")
nb_interFrames <- 25
nb_endDays <- 20
infected_per_point <- 1000
america_shift <- 40 # Shift America eastwards to reduce width of map
use_recovered <- FALSE
use_make_valid <- TRUE
outdir <- "/tmp/frames_corovir"
logfile <- "corovir_output.txt"
webserver_path <- "nyk:/var/www/nf/"
videofile <- sprintf("corovideo_%s.mp4", gsub("-", "", Sys.Date()))

sfc_shift <- function(geometry, x=0, y=0) {
	polygons <- sf::st_cast(geometry, "POLYGON")
	matrixList <- list()
	for (i in 1:nrow(polygons)) {
		tm <- matrix(rep(c(x, y), nrow(polygons[i,]$geometry[[1]][[1]])), ncol=2, byrow=TRUE)		
		matrixList[[i]] <- polygons[i,]$geometry[[1]][[1]] + tm
	}
	mpoly <- sf::st_multipolygon(list(matrixList))
	sf_mpoly <- sf::st_sf(sf::st_sfc(mpoly, crs=sf::st_crs(geometry)))
	if (!use_make_valid) return(sf_mpoly)
	return(sf::st_make_valid(sf_mpoly))
}

# Update local copy
git2r::pull("COVID-19")
# Load COVID data
covid <- read.csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv", stringsAsFactors=FALSE)
covid_rec <- read.csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_recovered_global.csv", stringsAsFactors=FALSE)
covid_dead <- read.csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv", stringsAsFactors=FALSE)
covdates <- as.Date(strptime(colnames(covid)[5:ncol(covid)], "X%m.%d.%y"))
# Check of there is new data
if (file.exists(".latest_covid_date.txt")) {
	ldate <- as.Date(scan(file=".latest_covid_date.txt", what="character"))
	if (max(covdates) <= ldate) writeLines("No new data!")
}
write(as.character(max(covdates)), file=".latest_covid_date.txt")

# Add "fake dates" at the end to not just stops at the last day of data
for (i in 1:nb_endDays) {
	covid <- cbind(covid, covid[,ncol(covid)])
	covid_rec	<- cbind(covid_rec, covid_rec[,ncol(covid_rec)])
	covid_dead	<- cbind(covid_dead, covid_dead[,ncol(covid_dead)])
}
colnames(covid)[(ncol(covid)-nb_endDays+1):ncol(covid)] <- strftime(seq.Date(from=max(covdates)+1, to=max(covdates)+nb_endDays, by=1), "X%m.%d.%y")
colnames(covid_rec)[(ncol(covid_rec)-nb_endDays+1):ncol(covid_rec)] <- strftime(seq.Date(from=max(covdates)+1, to=max(covdates)+nb_endDays, by=1), "X%m.%d.%y")
colnames(covid_dead)[(ncol(covid_dead)-nb_endDays+1):ncol(covid_dead)] <- strftime(seq.Date(from=max(covdates)+1, to=max(covdates)+nb_endDays, by=1), "X%m.%d.%y")
# Convert dates
covdates2 <- as.Date(strptime(colnames(covid)[5:ncol(covid)], "X%m.%d.%y"))
colnames(covid)[5:ncol(covid)] <- as.character(covdates2)
colnames(covid_rec)[5:ncol(covid_rec)] <- as.character(covdates2)
colnames(covid_dead)[5:ncol(covid_dead)] <- as.character(covdates2)

country_recode <- function(covid, world) {
	covid$country_name <- covid$Country.Region
	covid[covid$Province.State %in% world$name, "country_name"] <- covid[covid$Province.State %in% world$name, "Province.State"] 
	covid[covid$country_name == "US", "country_name"] <- "United States of America"
	covid[covid$country_name == "Korea, South", "country_name"] <- "South Korea"
	covid[covid$country_name == "Congo (Brazzaville)", "country_name"] <- "Congo"
	covid[covid$country_name == "Congo (Kinshasa)", "country_name"] <- "Congo"
	covid[covid$country_name == "Taiwan*", "country_name"] <- "Taiwan"
	#print(covid$country_name[!covid$country_name %in% world$name])
	covid[covid$country_name %in% world$name,]
}

# Prepare world map
world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")

# Shift America east
for (i in 1:nrow(world)) if (world[i, "continent"]$continent %in% c("North America", "South America")) world[i, "geometry"] <- sfc_shift(world[i, "geometry"], x=america_shift)

covid <- country_recode(covid, world)
covid_rec <- country_recode(covid_rec, world)
covid_dead <- country_recode(covid_dead, world)

country_points <- function(lat, long, country_name, nb_points) {
	if (world[world$name == country_name, "continent"]$continent %in% c("North America", "South America")) long <- long + america_shift
	polygons <- sf::st_cast(world[world$name == country_name, "geometry"], "POLYGON")
	polymatch <- as.vector(sf::st_within(sf::st_point(c(long, lat)), polygons, sparse = FALSE))
	if (sum(polymatch) == 0) {
		write(sprintf("%s: Coordiante not found in any geometry.", country_name), file=logfile, append=TRUE)
		return(NULL)
	}
	geom <- polygons[which(polymatch),]
	bbox <- sf::st_bbox(geom)
	nb_points2 <- nb_points * 5
	points <- data.frame(
		country=country_name,
		country_polygon=which(polymatch),
		long=rnorm(nb_points2, mean=long, sd=(bbox$xmax - bbox$xmin) / 5),
		lat=rnorm(nb_points2, mean=lat, sd=(bbox$ymax - bbox$ymin) / 5), 
		xvec=rnorm(nb_points2, mean=0, sd=0.15), 
		yvec=rnorm(nb_points2, mean=0, sd=0.15),
		status="infected",
		inside=TRUE,
		stringsAsFactors=FALSE
	)
	points_sf <- sf::st_as_sf(points, coords=c("long", "lat"), crs=sf::st_crs(geom))
	points2 <- points[sf::st_within(points_sf, geom, sparse=FALSE),]
	points2[1:nb_points,]
}

# Prepare map plot
worldPlot <- geom_sf(data=world, fill="black") 
viewBox <- coord_sf(xlim=c(-125+america_shift, 145), ylim=c(-35, 70), expand=FALSE)
pointColors <- scale_color_manual(values=c("green", "red"))
plotTheme <- theme(
	legend.position = "none", 
	axis.title.x=element_blank(), 
	axis.title.y=element_blank(), 
	plot.background=element_rect(fill = "black"), 
	panel.background = element_rect(fill = "midnightblue"), 
	plot.title=element_text(color = "lightgray")
)
p0 <- ggplot() + worldPlot + viewBox + plotTheme + pointColors

# Prepare output directory
if (file.exists(logfile)) unlink(logfile)
if (!dir.exists(outdir)) dir.create(outdir)
for (i in list.files(outdir, full.names=TRUE)) unlink(i)

covid$change <- 0
covid$change_dead <- 0
framenum <- 1
if (exists("cpt")) rm("cpt")
for (nowi in 1:length(covdates2)) {
	now <- covdates2[nowi]
	for (i in 1:nrow(covid)) {
		# Newly infected
		newinf <- covid[i, as.character(now)]
		if (is.na(newinf)) newinf <- 0
		if (as.character(now-1) %in% colnames(covid)) oldinf <- covid[i, as.character(now-1)] else oldinf <- 0
		newinf <- newinf - oldinf
		# Newly recovered
		newrec <- covid_rec[i, as.character(now)]
		if (is.na(newrec)) newrec <- 0
		if (as.character(now-1) %in% colnames(covid_rec)) oldrec <- covid_rec[i, as.character(now-1)] else oldrec <- 0
		if (is.na(oldrec)) oldrec <- 0
		newrec <- newrec - oldrec
		if (use_recovered) covid[i, "change"] <- covid[i, "change"] + newinf - newrec else covid[i, "change"] <- covid[i, "change"] + newinf
		if (covid[i, "change"] > 0) { # Add points to table
			nb_newpoints <- covid[i, "change"] %/% infected_per_point
			covid[i, "change"] <- covid[i, "change"] %% infected_per_point
			if (nb_newpoints > 0) {
				cp <- country_points(covid[i, "Lat"], covid[i, "Long"], covid[i, "country_name"], nb_newpoints)
				if (exists("cpt")) cpt <- rbind(cpt, cp) else cpt <- cp
				write(sprintf("%s: Infected %d * %d in %s (%d left)", now, nb_newpoints, infected_per_point, covid[i, "country_name"], covid[i, "change"]), file=logfile, append=TRUE)
			}
		}
		if (covid[i, "change"] < 0 & exists("cpt")) { # Remove points from table
			nb_recpoints <- (-1 * covid[i, "change"]) %/% infected_per_point
			covid[i, "change"] <- -1 * (-1 * covid[i, "change"]) %% infected_per_point
			if (nb_recpoints > 0) {
				coi <- cpt[cpt$country == covid[i, "country_name"] & cpt$status == "infected",]
				remrows <- rownames(coi)[1:nb_recpoints]
				cpt <- cpt[!rownames(cpt) %in% remrows,]
				write(sprintf("%s: Recovered %d * %d in %s (%d left)", now, nb_recpoints, infected_per_point, covid[i, "country_name"], covid[i, "change"]), file=logfile, append=TRUE)
			}
		}
		# Newly died
		newdead <- covid_dead[i, as.character(now)]
		if (is.na(newdead)) newdead <- 0
		if (as.character(now-1) %in% colnames(covid_dead)) olddead <- covid_dead[i, as.character(now-1)] else olddead <- 0
		newdead <- newdead - olddead
		covid[i, "change_dead"] <- covid[i, "change_dead"] + newdead
		if (covid[i, "change_dead"] > 0 & exists("cpt")) {
			nb_deadpoints <- covid[i, "change_dead"] %/% infected_per_point
			covid[i, "change_dead"] <- covid[i, "change_dead"] %% infected_per_point
			if (nb_deadpoints > 0) {  # Mark points from table as dead
				coi <- cpt[cpt$country == covid[i, "country_name"] & cpt$status == "infected",]
				deadrows <- rownames(coi)[1:nb_deadpoints]
				cpt[rownames(cpt) %in% deadrows, "status"] <- "dead"
				cpt[rownames(cpt) %in% deadrows, "xvec"] <- 0
				cpt[rownames(cpt) %in% deadrows, "yvec"] <- 0
				write(sprintf("%s: Died %d * %d in %s (%d left)", now, nb_deadpoints, infected_per_point, covid[i, "country_name"], covid[i, "change_dead"]), file=logfile, append=TRUE)
			}
		}
	}
	if (!exists("cpt")) next
	if (nrow(cpt) == 0) next
	write(sprintf("%s: Missing: %d of %d (%2.2f%%)", now, sum(is.na(cpt$lat)), nrow(cpt), 100 * sum(is.na(cpt$lat)) / nrow(cpt)), file=logfile, append=TRUE)
	print(table(cpt$country))
	cpt <- cpt[!is.na(cpt$lat),]
	cpt$status <- factor(cpt$status, levels=c("infected", "dead"))
	# Movement of points
	for (i in 1:nb_interFrames) {
		# Write map as PNG
		ofn <- sprintf("%s/frame%05d.png", outdir, framenum)
		framenum <- framenum + 1
		# Compute point size per country
		countryCount <- as.data.frame.table(table(cpt$country))
		countryCount$psize <- 4 - log(countryCount$Freq)
		countryCount[countryCount$psize < 1, "psize"] <- 1
		colnames(countryCount)[1] <- "country"
		cpt_psize <- plyr::join(cpt, countryCount, by="country")
		cptr <- cpt_psize[nrow(cpt_psize):1,] # To make the earliest infected (and dead) visible on top
		if (now > max(covdates)) nowk <- max(covdates) else nowk <- now
		p <- p0 + geom_point(aes(x=long, y=lat, color=status, size=psize), data=cptr, alpha=0.5) + ggtitle(nowk)
		png(filename=ofn, width=1600, height=750, bg="black")
		print(p)
		dev.off()
		# Move according to vector
		cpt$lat <- cpt$lat + cpt$xvec
		cpt$long <- cpt$long + cpt$yvec
		# Check if still inside country polygon
		for (country_name in unique(cpt$country)) {
			cpti <- cpt[cpt$country == country_name,]
			polygons <- sf::st_cast(world[world$name == country_name, "geometry"], "POLYGON")
			geom <- polygons[cpti$country_polygon[1],]
			points_sf <- sf::st_as_sf(cpti, coords=c("long", "lat"), crs=sf::st_crs(geom))
			pointmatch <- sf::st_within(points_sf, geom, sparse=FALSE)
			cpt[cpt$country == country_name, "inside"] <- pointmatch
		}
		# Move back inside polygon if outside
		cpt[cpt$inside == FALSE, "lat"] <- cpt[cpt$inside == FALSE, "lat"] - cpt[cpt$inside == FALSE, "xvec"]
		cpt[cpt$inside == FALSE, "long"] <- cpt[cpt$inside == FALSE, "long"] - cpt[cpt$inside == FALSE, "yvec"]
		# Randomly invert movement vector
		cpt[cpt$inside == FALSE, "xvec"] <- cpt[cpt$inside == FALSE, "xvec"] * (2 * round(runif(sum(!cpt$inside), min=0, max=1)) - 1)
		cpt[cpt$inside == FALSE, "yvec"] <- cpt[cpt$inside == FALSE, "yvec"] * (2 * round(runif(sum(!cpt$inside), min=0, max=1)) - 1)
		# Collision detection
		pDistMat <- sp::spDists(as.matrix(cpt[,c("long", "lat")])) # Create distance matrix
		touchPairs <- which(pDistMat < 0.25, arr.ind=TRUE) # Distance below approximately 25km (50km)
		touchIndex <- touchPairs[touchPairs[,1] != touchPairs[,2], 1] # Remove self distance
		# Randomly invert movement vector in case of collision
		if (length(touchIndex) > 0) {
			write(sprintf("%s: %d collisions", now, length(touchIndex)), file=logfile, append=TRUE)
			cpt[touchIndex, "xvec"] <- cpt[touchIndex, "xvec"] * (2 * round(runif(length(touchIndex), min=0, max=1)) - 1)
			cpt[touchIndex, "yvec"] <- cpt[touchIndex, "yvec"] * (2 * round(runif(length(touchIndex), min=0, max=1)) - 1)
		}
	}
}

# Create video from frames
cmd <- sprintf("ffmpeg -y -i %s/frame%%05d.png -c:v libx264 -strict -2 -pix_fmt yuv420p -f mp4 %s", outdir, videofile)
write(cmd, file=logfile, append=TRUE)
system(cmd)

# Upload to webserver
cmd <- sprintf("rsync -vrpe ssh %s %s", videofile, webserver_path)
write(cmd, file=logfile, append=TRUE)
system(cmd)

# HTML5 video player
html <- sprintf('<html><head><title>COVID-19 pandemic visualisation by Niklaus Fankhauser</title></head><body style="background:black;color:white">
<h1>COVID-19 pandemic visualisation by Niklaus Fankhauser</h1>
<video controls autoplay loop><source src="%s" type="video/mp4"></video> 
<p>The idea is to represent the spread of the pandemic as ideal gases inside each country.
Dots represents 1000 infected (green) or dead (red).</p>
<p>Based on data from <i>2019 Novel Coronavirus COVID-19 (2019-nCoV) Data Repository by Johns Hopkins CSSE.</i></p>
<p>Updated on %s</p>
</body></html>', videofile, Sys.time())
write(html, file="corona.html")

# Upload to webserver
cmd <- sprintf("rsync -vrpe ssh corona.html %s", webserver_path)
write(cmd, file=logfile, append=TRUE)
system(cmd)
unlink("corona.html")
