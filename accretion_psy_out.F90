module accretion
  use variable_precision, only : wp
  use passive_fields, only : rho
  use mphys_switches, only : active_rain, i_am4, i_am5, i_cfl, i_cfr, i_nl, i_ql, i_qr, isol, l_2mc, l_aacc, l_kk00, l_preventsmall, l_prf_cfrac, l_process
  use mphys_constants, only : fixed_cloud_number
  use mphys_parameters, only : hydro_params
  use process_routines, only : i_aacw, i_pracw, process_rate
  use thresholds, only : cfliq_small, ql_small, qr_small
  use sweepout_rate, only : sweepout
  use distributions, only : dist_lambda, dist_mu, dist_n0
  implicit none
  character(len = *), parameter, private :: ModuleName = 'ACCRETION'
  private

  public :: racw

  contains
  subroutine racw(ixy_inner, dt, qfields, cffields, aerofields, procs, params, aerosol_procs)
    use yomhook, only : dr_hook, lhook
    use parkind1, only : jpim, jprb
    integer(kind=jpim), parameter :: zhook_in = 0
    integer(kind=jpim), parameter :: zhook_out = 1
    integer, intent(in) :: ixy_inner
    real(kind=wp), intent(in) :: dt
    real(kind=wp), dimension(:,:), intent(in) :: qfields
    real(kind=wp), dimension(:,:), intent(in) :: cffields
    real(kind=wp), dimension(:,:), intent(in) :: aerofields
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:, :)
    TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:, :)
    type(hydro_params), intent(in) :: params
    real(kind=wp) :: dmass
    real(kind=wp) :: dnumber
    real(kind=wp) :: damass
    real(kind=wp) :: cloud_mass
    real(kind=wp) :: cloud_number
    real(kind=wp) :: rain_mass
    real(kind=wp) :: cf_liquid
    real(kind=wp) :: cf_rain
    real(kind=wp) :: mu
    real(kind=wp) :: n0
    real(kind=wp) :: lam
    logical, save :: l_kk_acw = .true.
    integer :: k
    CHARACTER(LEN = *), PARAMETER :: RoutineName = 'RACW'
    real(kind=jprb) :: zhook_handle
    integer :: loop_stop
    integer :: loop_stop_1

    !$acc routine vector
    loop_stop = UBOUND(qfields, 1)
    !$acc loop vector
    do k = 1, loop_stop, 1
      if (l_prf_cfrac) then
        if (cffields(k,i_cfl) > cfliq_small) then
          cf_liquid = cffields(k,i_cfl)
        else
          cf_liquid = cfliq_small
        end if
        if (cffields(k,i_cfr) > cfliq_small) then
          cf_rain = cffields(k,i_cfr)
        else
          cf_rain = cfliq_small
        end if
      else
        cf_liquid = 1.0_wp
        cf_rain = 1.0_wp
      end if
      cloud_mass = qfields(k,i_ql) / cf_liquid
      rain_mass = qfields(k,i_qr) / cf_rain
      if (l_2mc) then
        cloud_number = qfields(k,i_nl) / cf_liquid
      else
        cloud_number = fixed_cloud_number / cf_liquid
      end if
      if (cloud_mass * cf_liquid > ql_small .AND. rain_mass * cf_rain > qr_small) then
        if (l_kk_acw) then
          if (l_kk00) then
            dmass = MIN(0.9 * cloud_mass, 67.0 * (cloud_mass * rain_mass) ** 1.15)
          else
            dmass = MIN(0.9 * cloud_mass, 8.53 * cloud_mass ** 1.05 * rain_mass ** 0.98)
          end if
        else
          n0 = dist_n0(k,params%id,ixy_inner)
          mu = dist_mu(k,params%id,ixy_inner)
          lam = dist_lambda(k,params%id,ixy_inner)
          dmass = sweepout(n0,lam,mu,params,rho(k,ixy_inner)) * cloud_mass
        end if
        if (l_preventsmall .AND. dmass < qr_small) then
          dmass = 0.0
        end if
        if (l_2mc) then
          dnumber = dmass / (cloud_mass / cloud_number)
        end if
        if (l_prf_cfrac) then
          dmass = dmass * MIN(cf_liquid, cf_rain)
          dnumber = dnumber * MIN(cf_liquid, cf_rain)
          cloud_mass = cloud_mass * cf_liquid
        end if
        procs(i_ql,i_pracw%id)%column_data(k) = -dmass
        procs(i_qr,i_pracw%id)%column_data(k) = dmass
        if (l_2mc) then
          procs(i_nl,i_pracw%id)%column_data(k) = -dnumber
        end if
      end if
    enddo
    if (l_aacc .AND. l_process) then
      if (active_rain(isol)) then
        loop_stop_1 = UBOUND(qfields, 1)
        !$acc loop vector
        do k = 1, loop_stop_1, 1
          dmass = procs(i_qr,i_pracw%id)%column_data(k)
          damass = dmass / cloud_mass * aerofields(k,i_am4)
          aerosol_procs(i_am4,i_aacw%id)%column_data(k) = -damass
          aerosol_procs(i_am5,i_aacw%id)%column_data(k) = damass
        enddo
      end if
    end if

  end subroutine racw

end module accretion
