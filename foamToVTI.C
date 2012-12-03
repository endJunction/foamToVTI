
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

vtkDoubleArray*
createDoubleArray(const std::string& name, const size_t length, const size_t components = 1)
{
    vtkDoubleArray* a = vtkDoubleArray::New();
    a->SetName(name.c_str());
    a->SetNumberOfComponents(components);
    a->SetNumberOfTuples(length);
    a->Allocate(length);

    double* value = new double[components];
    for (size_t i = 0; i < components; i++)
        value[i] = 0;

    for (size_t i = 0; i < length; i++)
        a->InsertTuple(i, value);

    delete[] value;
    return a;
}

vtkUnsignedCharArray*
createBlankingArray(const size_t length)
{
    vtkUnsignedCharArray* a = vtkUnsignedCharArray::New();
    a->SetName("blanking");
    a->SetNumberOfComponents(1);
    a->SetNumberOfTuples(length);
    a->Allocate(length);

    for (size_t i = 0; i < length; i++)
        a->InsertTuple1(i, -1);

    return a;
}

vtkUniformGrid*
createUniformGrid(const boundBox& bounds, const double h,
    const bool createPointData = true)
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

    const vtkIdType nCells = image->GetNumberOfCells();
    image->GetCellData()->AddArray(createDoubleArray("U", nCells, 3));
    image->GetCellData()->GetVectors("U")->Delete();
    image->GetCellData()->AddArray(createDoubleArray("p", nCells, 1));
    image->GetCellData()->GetScalars("p")->Delete();
    image->GetCellData()->AddArray(createBlankingArray(nCells));
    image->GetCellData()->GetScalars("blanking")->Delete();

    image->GetCellData()->SetActiveVectors("U");
    image->GetCellData()->SetActiveScalars("p");

    if (!createPointData)
        return image;

    const vtkIdType nPoints = image->GetNumberOfPoints();
    image->GetPointData()->AddArray(createDoubleArray("U", nPoints, 3));
    image->GetPointData()->GetVectors("U")->Delete();
    image->GetPointData()->AddArray(createDoubleArray("p", nPoints, 1));
    image->GetPointData()->GetScalars("p")->Delete();
    image->GetPointData()->AddArray(createBlankingArray(nPoints));
    image->GetPointData()->GetScalars("blanking")->Delete();

    image->GetPointData()->SetActiveVectors("U");
    image->GetPointData()->SetActiveScalars("p");

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
    argList::addBoolOption
    (
        "internalPoints",
        "save point data too. Default is to save cell data only"
    );

    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"

    runTime.setTime(timeDirs.last(), timeDirs.size()-1);

    const bool saveInternalPoints = args.optionFound("internalPoints");

    vtkUniformGrid* image
        = createUniformGrid(mesh.bounds(), 1./32, saveInternalPoints);
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

        vtkDataArray* cU = image->GetCellData()->GetVectors("U");
        if (!cU)
        {
            std::cerr << "Cells' vector U not defined.\n";
            return EXIT_FAILURE;
        }

        vtkDataArray* cP = image->GetCellData()->GetScalars("p");
        if (!cP)
        {
            std::cerr << "Cells' scalars p not defined.\n";
            return EXIT_FAILURE;
        }


        for (label i = 0; i < mesh.nCells(); i++)
        {
            double coords[3]
                = { cellCenters[i].x(),
                    cellCenters[i].y(),
                    cellCenters[i].z() };
            int subId = 0;
            double pcoords[3];
            double weights[8];  // For 3D image cell has 8 points.

            const vtkIdType cellId
                = image->FindCell(coords, 0, 0, 1e-8, subId, pcoords, weights);

            cU->SetTuple(cellId, &(U[i][0]));
            cP->SetTuple(cellId, &(p[i]));
            image->UnBlankCell(cellId);
        }

        if (saveInternalPoints)
        {
            vtkDataArray* pU = image->GetPointData()->GetVectors("U");
            if (!pU)
            {
                std::cerr << "Points' vector U not defined.\n";
                return EXIT_FAILURE;
            }

            vtkDataArray* pP = image->GetPointData()->GetScalars("p");
            if (!pP)
            {
                std::cerr << "Points' scalars p not defined.\n";
                return EXIT_FAILURE;
            }

            for (label i = 0; i < mesh.nPoints(); i++)
            {
                const vtkIdType id = image->FindPoint(
                        pointCoords[i].x(), pointCoords[i].y(), pointCoords[i].z());

                pU->SetTuple(id, &(pointU[i][0]));
                pP->SetTuple(id, &pointP[i]);
                image->UnBlankPoint(id);
            }
        }
        writeImage(image, runTime.timeName());
    }

    Info<< endl;

    image->Delete();

    return EXIT_SUCCESS;
}


// ************************************************************************* //
