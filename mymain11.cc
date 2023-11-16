#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"
using namespace Pythia8;

int main() {
    // Generator
    Pythia pythia;

    // read events from POWHEG lhe output
    pythia.readString("Beams:frameType = 4");
    pythia.readString("Beams:LHEF = pwgevents.lhe");

    // Event settings
    pythia.readString("Main:numberOfEvents = 0");
    int nError = 10;
    pythia.readString("Main:timesAllowErrors = 10");

    // POWHEG settings
    pythia.readString("POWHEG:veto = 1");
    pythia.readString("POWHEG:pTdef = 1");
    pythia.readString("POWHEG:emitted = 0");
    pythia.readString("POWHEG:pTemt = 0");
    pythia.readString("POWHEG:pThard = 0");
    pythia.readString("POWHEG:vetoCount = 3");
    pythia.readString("POWHEG:MPIveto = 0");

    // Add in user hooks for shower vetoing.
    shared_ptr<PowhegHooks> powhegHooks;
    powhegHooks = make_shared<PowhegHooks>();
    pythia.setUserHooksPtr((UserHooksPtr)powhegHooks);

    // Initialise and list settings
    pythia.init();

    // Begin event loop; generate until end of LHEF file
    int iEvent = 0, iError = 0;
    while (true) {
        // Generate the next event
        if (!pythia.next()) {
            // If failure because reached end of file then exit event loop
            if (pythia.info.atEndOfFile()) break;

            // Otherwise count event failure and continue/exit as necessary
            cout << "Warning: event " << iEvent << " failed" << endl;
            if (++iError == nError) {
                cout << "Error: too many event failures... exiting" << endl;
                break;
            }

            continue;
        }

        // process the event

        ++iEvent;
    }

    pythia.stat();

    return 0;
}