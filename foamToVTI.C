
#include <vtkUniformGrid.h>
#include <vtkXMLImageDataWriter.h>

#include "fvCFD.H"

vtkUniformGrid*
createUniformGrid(const boundBox& bounds, const double h)
{
    vtkUniformGrid* image = vtkUniformGrid::New();

    const vector l = bounds.max() - bounds.min();
    const unsigned int nPointsX = std::ceil(l.x()/h);
    const unsigned int nPointsY = std::ceil(l.y()/h);
    const unsigned int nPointsZ = std::ceil(l.z()/h);
    const double dX = l.x()/(nPointsX-1);
    const double dY = l.y()/(nPointsY-1);
    const double dZ = l.z()/(nPointsZ-1);


    image->SetExtent(0, nPointsX-1, 0, nPointsY-1, 0, nPointsZ-1);
    image->SetSpacing(dX, dY, dZ);
    image->SetOrigin(bounds.min().x(), bounds.min().y(), bounds.min().z());

    return image;
}

void
writeImage(vtkUniformGrid * const image, const std::string& timeName)
{
    vtkXMLImageDataWriter* w = vtkXMLImageDataWriter::New();
    w->SetFileName((timeName + ".vti").c_str());
    w->SetDataModeToBinary();
    w->SetInputData(image);
    w->Write();
    w->Delete();
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"

    runTime.setTime(timeDirs.last(), timeDirs.size()-1);

    vtkUniformGrid* image = createUniformGrid(mesh.bounds(), 1./32);

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject pheader
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (!pheader.headerOk())
        {
            Info << "Unable to read p" << endl;
            return EXIT_FAILURE;
        }

        Info << "    Reading p" << endl;
        volScalarField p(pheader, mesh);

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        if (!Uheader.headerOk())
        {
            Info << "Unable to read U" << endl;
            return EXIT_FAILURE;
        }

        Info << "    Reading U" << endl;
        volVectorField U(Uheader, mesh);

        writeImage(image, runTime.timeName());
    }

    Info<< endl;

    return EXIT_SUCCESS;
}


// ************************************************************************* //
