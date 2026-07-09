//- Set alpha between 0 and 1 only at the interface. Remove "dust".
{

    Switch hardBoundAlpha
    (
        mixture.lookupOrDefault<Switch>("alphaInterfaceOnly",false)
    );

    if(hardBoundAlpha)
    {
        //- Reconstruction is owned by the isoAdvection advector (native
        //  interfaceProperties has none). hardBoundAlpha is only reached in
        //  the isoAdvector branch, so the advector is always available here.
        const boolList& interfaceCells(advector.surf().interfaceCell());

        forAll(interfaceCells,celli)
        {
            if(!interfaceCells[celli])
            {
                alpha1[celli] = pos0(alpha1[celli]-0.5);
            }
        }
    }

}
