import numpy as np
import resonance
import sys


class TestParticle:
    def __init__(self):  # Attributes defined
        self.Resonant = False
        self.ResonanceType = 'n:n'
        self.Name = 'N/A'
        self.ResonanceCenter = 999
        self.ResonanceAmplitude = 999
        self.AverageSMA = 999  # Average SemiMajor axis
        self.AverageEccentricity = 999
        self.AverageInclination = 999
        self.Kozai = False
        self.SMAamplitude = 999
        self.SMACenter = 999

    ############################################      FUNCTIONS      #################################################
    ############################################   DATA DISSECTION   #################################################
    def DataDissection(self, typeOfData):

        TestParticleSample = sys.argv[1]  # User to choose a test sample using terminal

        with open(TestParticleSample + "bary.out") as f:  # Counting number of lines
            for line, l in enumerate(f):
                pass
        NumberOfLines = line

        # Taking the test point's data from the .out file sequentially
        TestParticleTime, SemiMajorAxis, Eccentricity, Inclination, Ω, ω, AngularPosition, Longitude = np.genfromtxt(
            str(TestParticleSample) + "bary.out", unpack=True)
        λ = (Ω + ω + AngularPosition) % 360  # The Lambda for test particles
        Pomega = (Ω + ω) % 360  # The longitude if pericenter in degrees

        # Flags "Specific ones"
        IsItResonant = False  # Is it in resonance?
        ResonanceAmplitude = 999  # The Resonance Amplitude
        ResonanceCenter = 999  # The Resonance Center
        ResonanceName = 999  # The Resonance name "Ration"
        IsItKozai = False  # Is it Kozai resonance?
        SMAAmplitude = 999  # SemiMajor amplitude
        SMACenter = 999  # SemiMajor center

        # Flags "General ones"
        IsIt = False  # Resonance / Kozai ?
        Amplitude = 999  # Phi / SMA
        Center = 999  # Phi / SMA
        Name = 999  # Name of the test particle

        # list of resonances to check: pp and qq for pp:qq resonance
        pp = [2, 3, 3, 4, 4, 5, 5, 5, 5, 6, 7, 7, 7, 7, 8, 8, 9, 9, 9, 10]
        qq = [1, 1, 2, 1, 3, 1, 2, 3, 4, 1, 1, 2, 3, 4, 1, 3, 1, 2, 4, 1]

        for jj in np.arange(0, len(pp)):  # First Loop
            ResSemiMajorAxis = 30.1 * (float(pp[jj]) / float(qq[jj])) ** (
                    2. / 3.)  # Kepler's Third Law to calculate semimajor axis of the resonance

            # Searching within 2 AUs from the resonance center
            if IsIt == 0 and (ResSemiMajorAxis + 2) > np.average(SemiMajorAxis) > (ResSemiMajorAxis - 2):
                phi = (float(pp[jj]) * λ - float(qq[jj]) * Longitude - (float(pp[jj]) - float(qq[jj])) * Pomega) % 360

                AngleRange = np.arange(0, 360, 15)  # Array of angles 15 degrees increment each step
                Window = int(0)
                Loop = 0

                if typeOfData == 0:
                    # Dividing the timeline to 10 separate windows Detecting resonance on smaller scales
                    WindowStep = int(NumberOfLines / 10)

                    IsItArray = np.zeros(int(len(
                        phi) / WindowStep))  # Array of 10 binary elements to check for resonance each step '10%' set to zero
                    CenterArray = np.zeros(int(len(
                        phi) / WindowStep))  # Array of 10 binary elements to check the res angle each step '10%' set to zero

                    while Window + WindowStep < len(phi):

                        # Average of the semi-major axis from Current Window -> Next Window
                        WindowAverage = np.average(SemiMajorAxis[Window:Window + WindowStep])
                        if (ResSemiMajorAxis + 2) > WindowAverage > (
                                ResSemiMajorAxis - 2):  # Within 2 AUs of Window Average
                            WindowPhi = phi[Window:Window + WindowStep]  # Phi of next window
                            AnglePresent = np.zeros(len(AngleRange)) + 1
                            for step in np.arange(0, len(
                                    AngleRange) - 1):  # find out where the res angle doesn't go for 15 degrees, proxy for AnglePresent
                                if len(WindowPhi[
                                           (WindowPhi > AngleRange[step]) * (WindowPhi < (AngleRange[step + 1]))]) == 0:
                                    AnglePresent[step] = 0
                            IsItArray[Loop] = np.average(AnglePresent) * 180.
                            CenterArray[Loop] = np.average(
                                AnglePresent[AnglePresent != 0] * AngleRange[AnglePresent != 0])
                        else:
                            IsItArray[Loop] = 180.
                        Window += WindowStep  # Increment Window
                        Loop += 1  # Increment Loop
                    if len(IsItArray[
                               IsItArray < 180.]) > 8:  # If 8 out of 10 Windows classified as Resonant
                        IsIt = True
                        Amplitude = np.average(IsItArray)
                        Center = np.average(CenterArray)
                        Name = str(pp[jj]) + ':' + str(qq[jj])
                    else:
                        Amplitude = 999.  # Ask about the reason of this one

                else:
                    # If checking for Kozai, we only want one window
                    WindowStep = int(NumberOfLines)
                    IsItArray = np.zeros(int(len(
                        SemiMajorAxis) / WindowStep))  # For Kozai we check SMA
                    CenterArray = np.zeros(int(len(
                        SemiMajorAxis) / WindowStep))

                    while Window + WindowStep < len(SemiMajorAxis):
                        WindowSMA = SemiMajorAxis[Window:Window + WindowStep]  # Phi of next window
                        AnglePresent = np.zeros(len(AngleRange)) + 1

                        for step in np.arange(0, len(
                                AngleRange) - 1):  # find out where the res angle doesn't go for 15 degrees, proxy for AnglePresent
                            if len(WindowSMA[
                                       (WindowSMA > AngleRange[step]) * (WindowSMA < (AngleRange[step + 1]))]) == 0:
                                AnglePresent[step] = 0
                        IsItArray[Loop] = np.average(AnglePresent) * 180.
                        CenterArray[Loop] = np.average(
                            AnglePresent[AnglePresent != 0] * AngleRange[AnglePresent != 0])

                        Window += WindowStep  # Increment Window
                        Loop += 1  # Increment Loop
                    if len(IsItArray[
                               IsItArray < 180.]) > 8:  # If 8 out of 10 Windows classified as Resonant
                        IsIt = True
                        Amplitude = np.average(IsItArray)
                        Center = np.average(CenterArray)
                        Name = str(pp[jj]) + ':' + str(qq[jj])
                        break
                    else:
                        Amplitude = 999.  # Ask about the reason of this one

        if typeOfData == 0:
            IsItResonant = IsIt
            ResonanceAmplitude = Amplitude
            ResonanceCenter = Center
            ResonanceName = Name

            # Check if Resonance center is having large range
            MaxCenter = max(CenterArray)
            MinCenter = min(CenterArray)
            if (MaxCenter - MinCenter) > 100:
                IsItResonant = False

            self.Resonant = IsItResonant
            self.ResonanceAmplitude = ResonanceAmplitude
            self.ResonanceCenter = ResonanceCenter
            self.ResonanceType = ResonanceName

        else:
            IsItKozai = IsIt
            SMAAmplitude = Amplitude
            SMACenter = Center

            self.Kozai = IsItKozai
            self.SMAamplitude = SMAAmplitude
            self.SMACenter = SMACenter

        self.Name = TestParticleSample
        self.AverageEccentricity = np.average(Eccentricity)
        self.AverageInclination = np.average(Inclination)
        self.AverageSMA = np.average(SemiMajorAxis)

        return

    ############################################   IDENTIFY RESONANCE   ##############################################
    def IdentifyResonance(self):
        type = 0  # Indicated that the variable Resonant is what we want from DataDissection function
        self.DataDissection(type)
        if self.Resonant == True:
            type = 1  # Indicated that the variable Kozai is what we want from DataDissection function
            self.DataDissection(type)

    ##############################################      PRINT DATA      ##############################################
    def PrintData(self):  # To be changed to SetData at the end of the project
        print(self.Name)
        print("General Details:_____________________________________________________________")
        print(self.AverageSMA, self.AverageEccentricity, self.AverageInclination)
        print("Resonance Detail:____________________________________________________________")
        print(self.Name, self.AverageSMA, self.AverageEccentricity, self.AverageInclination,
              self.ResonanceType, self.ResonanceAmplitude, self.ResonanceCenter, self.Resonant)
        print("Kozai Detail:________________________________________________________________")
        print(self.Kozai, self.SMAamplitude, self.SMACenter)
