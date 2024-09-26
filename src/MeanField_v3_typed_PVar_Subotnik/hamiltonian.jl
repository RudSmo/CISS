using LinearAlgebra



function SubotnikModelElectronic(HamPar::Any)::Array
	Delta = HamPar.Delta
        A = HamPar.A
        B = HamPar.B
        return [Delta 0 A-1im*B 0;
		0 Delta 0 A+1im*B;
                A+1im*B 0 -Delta 0;
                0 A-1im*B 0 -Delta]

end

function SubotnikModel_Mode(HamPar::Any,i::Int64)::Array
	if i == 1
		X = 1.0
		Y = 0.0
	else
		X = 0.0
		Y = 0.0
	end
        Delta = HamPar.Delta
        A = HamPar.A
        B = HamPar.B
	
        return [Delta+X 0 A*X-1im*B*Y 0;
                0 Delta+X 0 A*X+1im*B*Y;
                A*X+1im*B*Y 0 -Delta-X 0;
                0 A*X-1im*B*Y 0 -Delta-X]
end

function SubotnikModel(HamPar::Any)::Array
	X = HamPar.X[1]
	Y = HamPar.X[2]
	Delta = HamPar.Delta
	A = HamPar.A
	B = HamPar.B
	return [Delta+X 0 A*X-1im*B*Y 0;
                0 Delta+X 0 A*X+1im*B*Y;
                A*X+1im*B*Y 0 -Delta-X 0;
                0 A*X-1im*B*Y 0 -Delta-X]
end

