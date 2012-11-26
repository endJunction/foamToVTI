
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkUniformGrid.h>
#include <vtkUnsignedCharArray.h>
#include <vtkXMLImageDataWriter.h>

#include "fvCFD.H"
#include "volPointInterpolation.H"

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

    vtkDoubleArray* imgU = vtkDoubleArray::New();
    imgU->SetName("U");
    imgU->SetNumberOfComponents(3);
    imgU->SetNumberOfTuples(image->GetNumberOfCells());
    imgU->Allocate(image->GetNumberOfCells());
    image->GetCellData()->SetVectors(imgU);

    vtkDoubleArray* imgP = vtkDoubleArray::New();
    imgP->SetName("p");
    imgP->SetNumberOfComponents(1);
    imgP->SetNumberOfTuples(image->GetNumberOfCells());
    imgP->Allocate(image->GetNumberOfCells());
    image->GetCellData()->SetScalars(imgP);

    vtkUnsignedCharArray* blankCells = vtkUnsignedCharArray::New();
    blankCells->SetName("cell blanking");
    blankCells->SetNumberOfComponents(1);
    blankCells->SetNumberOfTuples(image->GetNumberOfCells());
    blankCells->Allocate(image->GetNumberOfCells());
    image->GetCellData()->AddArray(blankCells);
    image->SetCellVisibilityArray(blankCells);

    if (image->GetCellBlanking())
        std::cout << "blanking set\n";
    else
        std::cout << "blanking NOT set\n";

    for (vtkIdType i = 0; i < image->GetNumberOfCells(); i++)
    {
        imgU->InsertTuple3(i, 0, 0, 0);
        imgP->InsertTuple1(i, 0);
        blankCells->InsertTuple1(i, -1); // All cells visible
    }

    vtkDoubleArray* pU = vtkDoubleArray::New();
    pU->SetName("U");
    pU->SetNumberOfComponents(3);
    pU->SetNumberOfTuples(image->GetNumberOfPoints());
    pU->Allocate(image->GetNumberOfPoints());
    image->GetPointData()->SetVectors(pU);

    vtkDoubleArray* pP = vtkDoubleArray::New();
    pP->SetName("p");
    pP->SetNumberOfComponents(1);
    pP->SetNumberOfTuples(image->GetNumberOfPoints());
    pP->Allocate(image->GetNumberOfPoints());
    image->GetPointData()->SetScalars(pP);

    vtkUnsignedCharArray* blankPoints = vtkUnsignedCharArray::New();
    blankPoints->SetName("blanking");
    blankPoints->SetNumberOfComponents(1);
    blankPoints->SetNumberOfTuples(image->GetNumberOfPoints());
    blankPoints->Allocate(image->GetNumberOfPoints());
    image->GetPointData()->AddArray(blankPoints);
    image->SetPointVisibilityArray(blankPoints);

    for (vtkIdType i = 0; i < image->GetNumberOfPoints(); i++)
    {
        pU->InsertTuple3(i, 0, 0, 0);
        pP->InsertTuple1(i, 0);
        blankPoints->InsertTuple1(i, -1);    // All points visible
    }
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
    std::cout << *image << std::endl;

    const volVectorField& cellCenters = mesh.C();
    const pointField& pointCoords = mesh.points();
    volPointInterpolation interpolateVolPoint(mesh);

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
        pointScalarField pointP = interpolateVolPoint.interpolate(p);

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
        pointVectorField pointU = interpolateVolPoint.interpolate(U);

        for (label cellI = 0; cellI < mesh.nCells(); cellI++)
        {
            double coords[3] = { cellCenters[cellI].x(), cellCenters[cellI].y(), cellCenters[cellI].z() };
            int subId = 0;
            double pcoords[3];
            double weights[8];  // For 3D image cell has 8 points.
            const vtkIdType cellId
                = image->FindCell(coords, 0, 0, 1e-8, subId, pcoords, weights);

            image->GetCellData()->GetVectors("U")->SetTuple(cellId, &(U[cellI][0]));
            image->GetCellData()->GetScalars("p")->SetTuple(cellId, &p[cellI]);
            image->UnBlankCell(cellId);
        }

        for (label i = 0; i < mesh.nPoints(); i++)
        {
            const vtkIdType id = image->FindPoint(
                    pointCoords[i].x(), pointCoords[i].y(), pointCoords[i].z());

            image->GetPointData()->GetVectors("U")->SetTuple(id, &(pointU[i][0]));
            image->GetPointData()->GetScalars("p")->SetTuple(id, &pointP[i]);
            image->UnBlankPoint(id);
        }
        writeImage(image, runTime.timeName());
    }

    Info<< endl;

    return EXIT_SUCCESS;
}


// ************************************************************************* //
