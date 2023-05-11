#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main() {
    Pythia pythia;
    Hist pT("top transverse momentum", 100, 0., 200.);
    Hist eta("top pseudorapidity", 100, -5., 5.);
    Hist mult("charged multiplicity", 100, -1, 399);

    pythia.readString("Top:gg2ttbar = on");
    pythia.readString("Top:qqbar2ttbar = on");
    pythia.readString("Next:numberShowEvent = 5");

    pythia.readString("Beams:eCM = 13000.");
    pythia.init();

    // event loop
    for (int iEvent = 0; iEvent < 1000; ++iEvent) {
        pythia.next();
        // iterate through event record
        int iTop = 0; // look for last top state before decay
        int multCount = 0;
        for (int i = 0; i < pythia.event.size(); ++i) {
            if (pythia.event[i].id() == 6) iTop = i;
            if (pythia.event[i].isFinal() && pythia.event[i].isCharged()) multCount++;
        }
        pT.fill( pythia.event[iTop].pT() );
        eta.fill( pythia.event[iTop].eta() );
        mult.fill( multCount );
    }

    pythia.stat();

    cout << pT << eta << mult;

    return 0;
}