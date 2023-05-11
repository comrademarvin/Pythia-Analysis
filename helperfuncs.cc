// Prints condensed version of particle mother-tree
void decayStatus(int particleIndex, Pythia8::Event eventObj, int decayLevel = 0) {
    if (particleIndex != 0) {  
        cout << decayLevel << " " << eventObj[particleIndex].name() << " (" << particleIndex << ", " << eventObj[particleIndex].id() << ", " << eventObj[particleIndex].status() << ")" << ":";  
        for (int momIndex: eventObj[particleIndex].motherList()) {
            cout << " " << eventObj[momIndex].name() << " (" << momIndex << ", " << eventObj[momIndex].id() << ", " << eventObj[momIndex].status() << ")";
        }
        cout << endl;

        decayLevel++;
        for (int momIndex: eventObj[particleIndex].motherList()) {
            if ((!eventObj[momIndex].isGluon()) && (abs(eventObj[momIndex].id())!=2212)) {
                decayStatus(momIndex, eventObj, decayLevel);
            }
        }
    }
}

// Prints condensed version of particle daughter-tree
void decayStatusDaughters(int particleIndex, Pythia8::Event eventObj, int decayLevel = 0) {
    if (particleIndex != 0) {  
        cout << decayLevel << " " << eventObj[particleIndex].name() << " (" << particleIndex << ", " << eventObj[particleIndex].id() << ", " << eventObj[particleIndex].status() << ")" << ":";  
        for (int daughterIndex: eventObj[particleIndex].daughterList()) {
            cout << " " << eventObj[daughterIndex].name() << " (" << daughterIndex << ", " << eventObj[daughterIndex].id() << ", " << eventObj[daughterIndex].status() << ")";
        }
        cout << endl;

        decayLevel++;
        for (int daughterIndex: eventObj[particleIndex].daughterList()) {
            if (!eventObj[daughterIndex].isFinal()) {
                decayStatusDaughters(daughterIndex, eventObj, decayLevel);
            }
        }
    }
}

// returns index of decay muon in event, returns -1 if no decay muon is found
int muonDecay(int index, Pythia8::Event eventObj) {
    if (abs(eventObj[index].id()) == 13) return index;

    for (int daughterIndex: eventObj[index].daughterList()) {
        int returnIndex = muonDecay(daughterIndex, eventObj);
        if (returnIndex != -1) return returnIndex;
    }

    return -1;
}