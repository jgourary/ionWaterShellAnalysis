package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"os"
	"strconv"
	"strings"
)

var shellDist = 3.6

var frameTime = 2e-12

var binWidth = 6e-12

func main() {
	filePath := "liquid-e100-v100.arc"
	outPath := "bins.txt"
	logPath := "log.txt"

	logFile, err := os.Create(logPath)
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + logPath)
		log.Fatal(err)
	}

	residences, frameCount := getResidences(filePath, logFile)
	bins := generateBins(frameCount)
	bins = doStatistics(residences, bins, logFile)
	writeHistogram(outPath, bins)
}

func writeHistogram(outPath string, bins []bin) {
	thisFile, err := os.Create(outPath)
	if err != nil {
		fmt.Println("Failed to create new fragment file: " + outPath)
		log.Fatal(err)
	}
	for _, bin := range bins {
		thisFile.WriteString(fmt.Sprintf("%e", bin.binBottom) + ", " + fmt.Sprintf("%e", bin.binCenter) + ", " + fmt.Sprintf("%e", bin.binTop) + ", " + strconv.Itoa(bin.count) + "\n")
	}
}

func generateBins(frameCount int) []bin {
	bins := make([]bin, frameCount)
	for i := 0; i < frameCount; i++ {
		var newBin bin
		newBin.binBottom = float64(i) * binWidth
		newBin.binCenter = newBin.binBottom + binWidth/2.0
		newBin.binTop = newBin.binBottom + binWidth
		bins[i] = newBin
	}
	return bins
}

func doStatistics(residences []residence, bins []bin, logFile *os.File) []bin {
	if len(residences) == 0 {
		fmt.Println("No completed residences found")
	} else {
		// get average
		avgResidenceTime := 0.0
		for _, residence := range residences {
			avgResidenceTime += residence.residenceTime
		}
		avgResidenceTime = avgResidenceTime / float64(len(residences))
		fmt.Println("Average residence time = " + fmt.Sprintf("%e", avgResidenceTime))

		// get histogram
		for i, residence := range residences {
			logFile.WriteString("Residence #" + strconv.Itoa(i) + ": t = " + fmt.Sprintf("%e", residence.residenceTime) + "(" + strconv.Itoa(residence.firstFrame) + "-" + strconv.Itoa(residence.lastFrame) + ")\n")
			for i := 0; i < len(bins); i++ {
				if residence.residenceTime < bins[i].binTop {
					bins[i].count++
					// fmt.Println("Found bin!")
					break
				}
			}
		}
	}
	binCount := 0
	for _, bin := range bins {
		binCount += bin.count
	}
	// fmt.Println("Bin count = " + strconv.Itoa(binCount))
	return bins
}

func getResidences(filePath string, logFile *os.File) ([]residence, int) {
	// open file
	file, err := os.Open(filePath)
	if err != nil {
		fmt.Println("Failed to open molecule file: " + filePath)
		log.Fatal(err)
	}

	scanner := bufio.NewScanner(file)
	atoms := make(map[int]*atom)

	var completedResidences []residence
	currentResidences := make(map[int]residence)

	structureCenter := []float64{0.0, 0.0, 0.0}
	var lastShellOxygens []int
	var currentShellOxygens []int
	frameCount := 0
	oxygenCount := 0

	for scanner.Scan() {
		line := scanner.Text()
		fields := strings.Fields(line)
		if len(fields) > 5 && !strings.Contains(line, "30.000000   30.000000   30.000000") {
			newAtom := line2atom(fields)
			if newAtom.id == 1 {
				structureCenter = []float64{newAtom.x, newAtom.y, newAtom.z}
			} else if newAtom.element == "O" && dist2center(structureCenter, newAtom) < shellDist {
				currentShellOxygens = append(currentShellOxygens, newAtom.id)
			}
			atoms[newAtom.id] = &newAtom
		} else if len(fields) == 1 {
			frameCount++

			// create new residence structures for entering atoms
			enteringAtoms := enteringAtoms(lastShellOxygens, currentShellOxygens)
			for _, atomID := range enteringAtoms {
				var newResidence residence
				newResidence.atomID = atomID
				newResidence.firstFrame = frameCount
				currentResidences[atomID] = newResidence
			}

			// remove and save residence structures for leaving atoms
			leavingAtoms := leavingAtoms(lastShellOxygens, currentShellOxygens)
			for _, atomID := range leavingAtoms {
				completedResidence := currentResidences[atomID]
				completedResidence.lastFrame = frameCount
				completedResidence.residenceTime = float64(completedResidence.lastFrame-completedResidence.firstFrame) * frameTime
				completedResidences = append(completedResidences, completedResidence)

				delete(currentResidences, atomID)

			}
			// print
			if len(leavingAtoms) > 0 || len(enteringAtoms) > 0 {
				logFile.WriteString("Frame number " + strconv.Itoa(frameCount) + ":\n")
				logFile.WriteString("Current Atoms:\n")
				logFile.WriteString(intArray2string(currentShellOxygens))
				if len(leavingAtoms) > 0 {
					logFile.WriteString("Atoms Leaving\n")
					logFile.WriteString(intArray2string(leavingAtoms))
				}
				if len(enteringAtoms) > 0 {
					logFile.WriteString("Atoms Entering\n")
					logFile.WriteString(intArray2string(enteringAtoms))
				}
			}

			// clear map
			atoms = make(map[int]*atom)
			oxygenCount += len(currentShellOxygens)
			lastShellOxygens = currentShellOxygens
			currentShellOxygens = make([]int, 0)

		}
	}
	oxygenFloat := float64(oxygenCount) / float64(frameCount)
	fmt.Println("Average Oxygens in First Shell = " + fmt.Sprintf("%e", oxygenFloat))
	return completedResidences, frameCount
}

func intArray2string(ints []int) string {
	string := ""
	for _, int := range ints {
		string += strconv.Itoa(int) + ", "
	}
	string += "\n"
	return string
}
func enteringAtoms(lastShell []int, currentShell []int) []int {
	var enteringAtoms []int
	for _, atomID := range currentShell {
		if !contains(lastShell, atomID) {
			enteringAtoms = append(enteringAtoms, atomID)
		}
	}
	return enteringAtoms
}

func leavingAtoms(lastShell []int, currentShell []int) []int {
	var leavingAtoms []int
	for _, atomID := range lastShell {
		if !contains(currentShell, atomID) {
			leavingAtoms = append(leavingAtoms, atomID)
		}
	}
	return leavingAtoms
}

func contains(ints []int, int int) bool {
	for _, v := range ints {
		if v == int {
			return true
		}
	}
	return false
}

func dist2center(c []float64, a atom) float64 {
	dist := math.Sqrt(math.Pow(c[0]-a.x, 2) + math.Pow(c[1]-a.y, 2) + math.Pow(c[2]-a.z, 2))
	return dist
}

func line2atom(fields []string) atom {
	var newAtom atom
	newAtom.id, _ = strconv.Atoi(fields[0])
	newAtom.element = fields[1]
	newAtom.x, _ = strconv.ParseFloat(fields[2], 64)
	newAtom.y, _ = strconv.ParseFloat(fields[3], 64)
	newAtom.z, _ = strconv.ParseFloat(fields[4], 64)
	newAtom.atomType = fields[5]
	return newAtom
}

type residence struct {
	atomID        int
	firstFrame    int
	lastFrame     int
	residenceTime float64
}

type bin struct {
	binTop    float64
	binBottom float64
	binCenter float64
	count     int
}

type atom struct {
	element  string
	x        float64
	y        float64
	z        float64
	atomType string
	id       int
}
