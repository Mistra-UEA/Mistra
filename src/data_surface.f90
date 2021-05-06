!
! Copyright 2015-2017 the Authors
!
! Licensed under the EUPL, Version 1.1 only (the "Licence");
!
! You may not use this work except in compliance with the Licence.
! You may obtain a copy of the Licence at:
!   https://joinup.ec.europa.eu/software/page/eupl
!
! Unless required by applicable law or agreed to in writing,
! software distributed under the Licence is distributed on an
! "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
! either express or implied.
!
! See the Licence for the specific language governing permissions
! and limitations under the Licence.


module data_surface

! Description :
! -----------
  ! Global declaration of model variables for surface and canopy module

! Author :
! ------
  ! Werner Schneider (original version)
  ! Josue Bock (rewritten for this version of Mistra)

! Declarations:
! ------------
! Modules used:
  use precision, only : &
       dp                    ! double precision kind

  implicit none

  public ! confirm this property with an explicit statement
  save   ! not mandatory if called from main program

! surface temperature for water surface
!--------------------------------------
  real (kind=dp) :: &
       tw           ! water surface temperature                          (K)

end module data_surface
